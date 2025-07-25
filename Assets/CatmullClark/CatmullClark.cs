using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using UnityEditor;
using Unity.Mathematics;

public static class CatmullClark
{
	// crease data
	public struct cc_Crease
	{
		public int nextID;
		public int prevID;
		public float sharpness;
	};

	// generic halfedge data
	public struct cc_Halfedge
	{
		public int twinID;
		public int nextID;
		public int prevID;
		public int faceID;
		public int edgeID;
		public int vertexID;
		public int uvID;
	};

	// specialized halfedge data for quad-only meshes
	public struct cc_Halfedge_Quad
	{
		public int twinID;
		public int edgeID;
		public int vertexID;
		public int uvID;
	};

	// mesh data-structure
	public class cc_Mesh
	{
		public int vertexCount;
		public int uvCount;
		public int halfedgeCount;
		public int edgeCount;
		public int faceCount;
		public int[] vertexToHalfedgeIDs;
		public int[] edgeToHalfedgeIDs;
		public int[] faceToHalfedgeIDs;
		public float3[] vertexPoints;
		public float2[] uvs;
		public cc_Halfedge[] halfedges;
		public cc_Crease[] creases;
	};

	// subd data-structure
	public class cc_Subd
	{
		public cc_Mesh cage;
		public float3[] vertexPoints;
		public cc_Halfedge_Quad[] halfedges;
		public cc_Crease[] creases;
		public int maxDepth;
	};


	/*******************************************************************************
	 * UV Encoding / Decoding routines
	 *
	 */
	public static float2 cc__DecodeUv(int uvEncoded)
	{
		uint tmp = (uint)uvEncoded;
		float2 uv = new float2(
			((tmp >> 0) & 0xFFFF) / 65535.0f,
			((tmp >> 16) & 0xFFFF) / 65535.0f
		);

		return uv;
	}

	public static int cc__EncodeUv(float2 uv)
	{
		uint u = (uint)math.round(uv[0] * 65535.0f);
		uint v = (uint)math.round(uv[1] * 65535.0f);
		uint tmp = ((u & 0xFFFFu) | ((v & 0xFFFFu) << 16));

		return (int)tmp;
	}


	/*******************************************************************************
	 * FaceCount -- Returns the number of faces
	 *
	 */
	public static int ccm_FaceCount(cc_Mesh mesh)
	{
		return mesh.faceCount;
	}


	/*******************************************************************************
	 * EdgeCount -- Returns the number of edges
	 *
	 */
	public static int ccm_EdgeCount(cc_Mesh mesh)
	{
		return mesh.edgeCount;
	}


	/*******************************************************************************
	 * CreaseCount -- Returns the number of creases
	 *
	 */
	public static int ccm_CreaseCount(cc_Mesh mesh)
	{
		return ccm_EdgeCount(mesh);
	}


	/*******************************************************************************
	 * HalfedgeCount -- Returns the number of halfedges
	 *
	 */
	public static int ccm_HalfedgeCount(cc_Mesh mesh)
	{
		return mesh.halfedgeCount;
	}


	/*******************************************************************************
	 * VertexCount -- Returns the number of vertices
	 *
	 */
	public static int ccm_VertexCount(cc_Mesh mesh)
	{
		return mesh.vertexCount;
	}


	/*******************************************************************************
	 * UvCount -- Returns the number of uvs
	 *
	 */
	public static int ccm_UvCount(cc_Mesh mesh)
	{
		return mesh.uvCount;
	}


	/*******************************************************************************
	 * FaceCountAtDepth -- Returns the number of faces at a given subdivision depth
	 *
	 * The number of faces follows the rule
	 *          F^{d+1} = H^d
	 * Therefore, the number of halfedges at a given subdivision depth d>= 0 is
	 *          F^d = 4^{d - 1} H^0,
	 * where H0 denotes the number of half-edges of the control cage.
	 *
	 */
	public static int ccm_FaceCountAtDepth_Fast(cc_Mesh cage, int depth)
	{
		Debug.Assert(depth > 0);
		int H0 = ccm_HalfedgeCount(cage);

		return (H0 << ((depth - 1) << 1));
	}

	public static int ccm_FaceCountAtDepth(cc_Mesh cage, int depth)
	{
		if (depth == 0)
		{
			return ccm_FaceCount(cage);
		}
		else
		{
			return ccm_FaceCountAtDepth_Fast(cage, depth);
		}
	}


	/*******************************************************************************
	 * EdgeCountAtDepth -- Returns the number of edges at a given subdivision depth
	 *
	 * The number of edges follows the rule
	 *          E^{d+1} = 2 E^d + H^d
	 * Therefore, the number of edges at a given subdivision depth d>= 0 is
	 *          E^d = 2^{d - 1} (2 E^0 + (2^d - 1) H^0),
	 * where H0 and E0 respectively denote the number of half-edges and edges
	 * of the control cage.
	 *
	 */
	public static int ccm_EdgeCountAtDepth_Fast(cc_Mesh cage, int depth)
	{
		Debug.Assert(depth > 0);
		int E0 = ccm_EdgeCount(cage);
		int H0 = ccm_HalfedgeCount(cage);
		int tmp = (int)(~(0xFFFFFFFF << depth)); // (2^d - 1) // TODO: CHECK UINT->INT CONVERSION

		return ((E0 << 1) + (tmp * H0)) << (depth - 1);
	}

	public static int ccm_EdgeCountAtDepth(cc_Mesh cage, int depth)
	{
		if (depth == 0)
		{
			return ccm_EdgeCount(cage);
		}
		else
		{
			return ccm_EdgeCountAtDepth_Fast(cage, depth);
		}
	}


	/*******************************************************************************
	 * HalfedgeCountAtDepth -- Returns the number of halfedges at a given subd depth
	 *
	 * The number of halfedges is multiplied by 4 at each subdivision step.
	 * Therefore, the number of halfedges at a given subdivision depth d>= 0 is
	 *          4^d H0,
	 * where H0 denotes the number of half-edges of the control cage.
	 *
	 */
	public static int ccm_HalfedgeCountAtDepth(cc_Mesh cage, int depth)
	{
		int H0 = ccm_HalfedgeCount(cage);

		return H0 << (depth << 1);
	}


	/*******************************************************************************
	 * CreaseCountAtDepth -- Returns the number of creases at a given subd depth
	 *
	 * The number of creases is multiplied by 2 at each subdivision step.
	 * Therefore, the number of halfedges at a given subdivision depth d>= 0 is
	 *          2^d C0,
	 * where C0 denotes the number of creases of the control cage.
	 *
	 */
	public static int ccm_CreaseCountAtDepth(cc_Mesh cage, int depth)
	{
		int C0 = ccm_CreaseCount(cage);

		return C0 << depth;
	}


	/*******************************************************************************
	 * VertexCountAtDepth -- Returns the number of vertex points at a given subd depth
	 *
	 * The number of vertices follows the rule
	 *          V^{d+1} = V^d + E^d + F^d
	 * For a quad mesh, the number of vertices at a given subdivision depth d>= 0 is
	 *          V^d = V0 + (2^{d} - 1) E0 + (2^{d} - 1)^2 F0,
	 * where:
	 * - V0 denotes the number of vertices of the control cage
	 * - E0 denotes the number of edges of the control cage
	 * - F0 denotes the number of faces of the control cage
	 * Note that since the input mesh may contain non-quad faces, we compute
	 * the first subdivision step by hand and then apply the formula.
	 *
	 */
	public static int ccm_VertexCountAtDepth_Fast(cc_Mesh cage, int depth)
	{
		Debug.Assert(depth > 0);
		int V0 = ccm_VertexCount(cage);
		int F0 = ccm_FaceCount(cage);
		int E0 = ccm_EdgeCount(cage);
		int H0 = ccm_HalfedgeCount(cage);
		int F1 = H0;
		int E1 = 2 * E0 + H0;
		int V1 = V0 + E0 + F0;
		int tmp = (int)(~(0xFFFFFFFF << (depth - 1))); // 2^{d-1} - 1 // TODO: CHECK UINT->INT CONVERSION

		return V1 + tmp * (E1 + tmp * F1);
	}

	public static int ccm_VertexCountAtDepth(cc_Mesh cage, int depth)
	{
		if (depth == 0)
		{
			return ccm_VertexCount(cage);
		}
		else
		{
			return ccm_VertexCountAtDepth_Fast(cage, depth);
		}
	}


	/*******************************************************************************
	 * UvCountAtDepth -- Returns the number of vertices uvs at a given subd depth
	 *
	 * The uvs are duplicated for each halfedge. We do this because uvs are
	 * discontinuous wrt to vertex points in the general case and so induce a
	 * different topology than that produced by the vertex points. Since uvs don't
	 * need much precision we encode them within the halfedge.uvID field using
	 * 16-bit precision for both components (see the cc__EncodeUV routine).
	 * So here we simply return the number of halfedges
	 *
	 */
	public static int ccm_UvCountAtDepth(cc_Mesh cage, int depth)
	{
		return ccm_HalfedgeCountAtDepth(cage, depth);
	}


	/*******************************************************************************
	 * Halfedge data accessors
	 *
	 */
	public static cc_Halfedge ccm__Halfedge(cc_Mesh mesh, int halfedgeID)
	{
		return mesh.halfedges[halfedgeID];
	}

	public static int ccm_HalfedgeTwinID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).twinID;
	}

	public static int ccm_HalfedgeNextID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).nextID;
	}

	public static int ccm_HalfedgePrevID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).prevID;
	}

	public static int ccm_HalfedgeVertexID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).vertexID;
	}

	public static int ccm_HalfedgeUvID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).uvID;
	}

	public static int ccm_HalfedgeEdgeID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).edgeID;
	}

	public static int ccm_HalfedgeFaceID(cc_Mesh mesh, int halfedgeID)
	{
		return ccm__Halfedge(mesh, halfedgeID).faceID;
	}

	public static float ccm_HalfedgeSharpness(cc_Mesh mesh, int halfedgeID)
	{
		return ccm_CreaseSharpness(mesh, ccm_HalfedgeEdgeID(mesh, halfedgeID));
	}

	public static float3 ccm_HalfedgeVertexPoint(cc_Mesh mesh, int halfedgeID)
	{
		return ccm_VertexPoint(mesh, ccm_HalfedgeVertexID(mesh, halfedgeID));
	}

	public static float2 ccm_HalfedgeVertexUv(cc_Mesh mesh, int halfedgeID)
	{
		return ccm_Uv(mesh, ccm_HalfedgeUvID(mesh, halfedgeID));
	}

	public static cc_Crease ccm__Crease(cc_Mesh mesh, int edgeID)
	{
		return mesh.creases[edgeID];
	}

	public static int ccm_CreaseNextID(cc_Mesh mesh, int edgeID)
	{
		return ccm__Crease(mesh, edgeID).nextID;
	}

	public static int ccm_CreasePrevID(cc_Mesh mesh, int edgeID)
	{
		return ccm__Crease(mesh, edgeID).prevID;
	}

	public static float ccm_CreaseSharpness(cc_Mesh mesh, int edgeID)
	{
		return ccm__Crease(mesh, edgeID).sharpness;
	}

	public static int ccm_HalfedgeFaceID_Quad(int halfedgeID)
	{
		return halfedgeID >> 2;
	}

	public static int ccm__ScrollFaceHalfedgeID_Quad(int halfedgeID, int direction)
	{
		int baseID = 3;
		int localID = (halfedgeID & baseID) + direction;

		return (halfedgeID & ~baseID) | (localID & baseID);
	}

	public static int ccm_HalfedgeNextID_Quad(int halfedgeID)
	{
		return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, +1);
	}

	public static int ccm_HalfedgePrevID_Quad(int halfedgeID)
	{
		return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, -1);
	}


	/*******************************************************************************
	 * Vertex data accessors
	 *
	 */
	public static float3 ccm_VertexPoint(cc_Mesh mesh, int vertexID)
	{
		return mesh.vertexPoints[vertexID];
	}

	public static float2 ccm_Uv(cc_Mesh mesh, int uvID)
	{
		return mesh.uvs[uvID];
	}


	/*******************************************************************************
	 * VertexToHalfedgeID -- Returns a halfedge ID that carries a given vertex
	 *
	 */
	public static int ccm_VertexPointToHalfedgeID(cc_Mesh mesh, int vertexID)
	{
		return mesh.vertexToHalfedgeIDs[vertexID];
	}


	/*******************************************************************************
	 * EdgeToHalfedgeID -- Returns a halfedge associated with a given edge
	 *
	 */
	public static int ccm_EdgeToHalfedgeID(cc_Mesh mesh, int edgeID)
	{
		return mesh.edgeToHalfedgeIDs[edgeID];
	}


	/*******************************************************************************
	 * FaceToHalfedgeID -- Returns a halfedge associated with a given face
	 *
	 */
	public static int ccm_FaceToHalfedgeID(cc_Mesh mesh, int faceID)
	{
		return mesh.faceToHalfedgeIDs[faceID];
	}

	public static int ccm_FaceToHalfedgeID_Quad(int faceID)
	{
		return faceID << 2;
	}


	/*******************************************************************************
	 * Vertex Halfedge Iteration
	 *
	 */
	public static int ccm_NextVertexHalfedgeID(cc_Mesh mesh, int halfedgeID)
	{
		int twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);

		return twinID >= 0 ? ccm_HalfedgeNextID(mesh, twinID) : -1;
	}

	public static int ccm_PrevVertexHalfedgeID(cc_Mesh mesh, int halfedgeID)
	{
		int prevID = ccm_HalfedgePrevID(mesh, halfedgeID);

		return ccm_HalfedgeTwinID(mesh, prevID);
	}


	/*******************************************************************************
	 * Create -- Allocates memory for a mesh of given vertex and halfedge count
	 *
	 */
	public static cc_Mesh ccm_Create(
			int vertexCount,
			int uvCount,
			int halfedgeCount,
			int edgeCount,
			int faceCount)
	{
		cc_Mesh mesh = new cc_Mesh();

		mesh.vertexCount = vertexCount;
		mesh.uvCount = uvCount;
		mesh.halfedgeCount = halfedgeCount;
		mesh.edgeCount = edgeCount;
		mesh.faceCount = faceCount;
		mesh.vertexToHalfedgeIDs = new int[vertexCount];
		mesh.edgeToHalfedgeIDs = new int[edgeCount];
		mesh.faceToHalfedgeIDs = new int[faceCount];
		mesh.halfedges = new cc_Halfedge[halfedgeCount];
		mesh.creases = new cc_Crease[edgeCount];
		mesh.vertexPoints = new float3[vertexCount];
		mesh.uvs = new float2[uvCount];

		return mesh;
	}


	/*******************************************************************************
	 * CreateFromGPU -- Initialize a cc_Mesh from a cc_MeshGPU from GPU library
	 *
	 */
	public static cc_Mesh ccm_CreateFromCPU(CatmullClarkGPU.cc_MeshGPU meshGPU)
	{
		cc_Mesh mesh = ccm_Create(meshGPU.vertexCount,
						  meshGPU.uvCount,
						  meshGPU.halfedgeCount,
						  meshGPU.edgeCount,
						  meshGPU.faceCount);

		meshGPU.vertexToHalfedgeIDs.GetData(mesh.vertexToHalfedgeIDs);
		meshGPU.edgeToHalfedgeIDs.GetData(mesh.edgeToHalfedgeIDs);
		meshGPU.faceToHalfedgeIDs.GetData(mesh.faceToHalfedgeIDs);
		meshGPU.halfedges.GetData(mesh.halfedges);
		meshGPU.creases.GetData(mesh.creases);
		meshGPU.vertexPoints.GetData(mesh.vertexPoints);
		meshGPU.uvs.GetData(mesh.uvs);

		return mesh;
	}


	/*******************************************************************************
	 * FaceCountAtDepth -- Returns the accumulated number of faces up to a given subdivision depth
	 *
	 */
	public static int ccs_CumulativeFaceCountAtDepth(cc_Mesh cage, int depth)
	{
		return ccs_CumulativeHalfedgeCountAtDepth(cage, depth) >> 2;
	}

	public static int ccs_CumulativeFaceCount(cc_Subd subd)
	{
		return ccs_CumulativeFaceCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * EdgeCountAtDepth -- Returns the accumulated number of edges up to a given subdivision depth
	 *
	 */
	public static int ccs_CumulativeEdgeCountAtDepth(cc_Mesh cage, int depth)
	{
		Debug.Assert(depth >= 0);
		int H0 = ccm_HalfedgeCount(cage);
		int E0 = ccm_EdgeCount(cage);
		int H1 = H0 << 2;
		int E1 = (E0 << 1) + H0;
		int D = depth;
		int A = (int)(~(0xFFFFFFFF << D)); //  2^{d} - 1 // TODO : CHECK UINT->INT CONVERSION

		return (A * (6 * E1 + A * H1 - H1)) / 6;
	}

	public static int ccs_CumulativeEdgeCount(cc_Subd subd)
	{
		return ccs_CumulativeEdgeCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * HalfedgeCount -- Returns the total number of halfedges stored by the subd
	 *
	 * The number of halfedges is multiplied by 4 at each subdivision step.
	 * It follows that the number of half-edges is given by the formula
	 *    H = H0 x sum_{d=0}^{D} 4^d
	 *      = H0 (4^{D+1} - 1) / 3
	 * where D denotes the maximum subdivision depth and H0 the number of
	 * halfedges in the control mesh.
	 *
	 */
	public static int ccs_CumulativeHalfedgeCountAtDepth(cc_Mesh cage, int maxDepth)
	{
		Debug.Assert(maxDepth >= 0);
		int D = maxDepth;
		int H0 = ccm_HalfedgeCount(cage);
		int H1 = H0 << 2;
		int tmp = (int)(~(0xFFFFFFFF << (D << 1))); // (4^D - 1) // TODO : CHECK UINT->INT CONVERSION

		return (H1 * tmp) / 3;
	}

	public static int ccs_CumulativeHalfedgeCount(cc_Subd subd)
	{
		return ccs_CumulativeHalfedgeCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * CreaseCount -- Returns the total number of creases stored by the subd
	 *
	 * The number of creases is multiplied by 2 at each subdivision step.
	 * It follows that the number of half-edges is given by the formula
	 *    C = C0 x sum_{d=0}^{D} 2^d
	 *      = C0 (2^{D+1} - 1)
	 * where D denotes the maximum subdivision depth and C0 the number of
	 * creases in the control mesh.
	 *
	 */
	public static int ccs_CumulativeCreaseCountAtDepth(cc_Mesh cage, int maxDepth)
	{
		Debug.Assert(maxDepth >= 0);
		int D = maxDepth;
		int C0 = ccm_CreaseCount(cage);
		int C1 = C0 << 1;
		int tmp = (int)~(0xFFFFFFFF << D); // (2^D - 1) // TODO : CHECK UINT->INT CONVERSION

		return (C1 * tmp);
	}

	public static int ccs_CumulativeCreaseCount(cc_Subd subd)
	{
		return ccs_CumulativeCreaseCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * CumulativeVertexCount -- Returns the total number of vertices computed by the subd
	 *
	 * The number of vertices increases according to the following formula at
	 * each subdivision step:
	 *  Vd+1 = Fd + Ed + Vd
	 * It follows that the number of vertices is given by the formula
	 *   Vd = F0 + (2^(d+1) - 1) E0 +
	 *      = 4 H0 (4^D - 1) / 3
	 * where D denotes the maximum subdivision depth and H0 the number of
	 * halfedges in the control mesh
	 *
	 */
	public static int ccs_CumulativeVertexCountAtDepth(cc_Mesh cage, int depth)
	{
		Debug.Assert(depth >= 0);
		int V0 = ccm_VertexCount(cage);
		int F0 = ccm_FaceCount(cage);
		int E0 = ccm_EdgeCount(cage);
		int H0 = ccm_HalfedgeCount(cage);
		int F1 = H0;
		int E1 = 2 * E0 + H0;
		int V1 = V0 + E0 + F0;
		int D = depth;
		int A = (int)(~(0xFFFFFFFF << (D)));     //  2^{d} - 1 // TODO : CHECK UINT->INT CONVERSION
		int B = (int)(~(0xFFFFFFFF << (D << 1)) / 3); // (4^{d} - 1) / 3 // TODO : CHECK UINT->INT CONVERSION

		return A * (E1 - (F1 << 1)) + B * F1 + D * (F1 - E1 + V1);
	}

	public static int ccs_CumulativeVertexCount(cc_Subd subd)
	{
		return ccs_CumulativeVertexCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * MaxDepth -- Retrieve the maximum subdivision depth of the subd
	 *
	 */
	public static int ccs_MaxDepth(cc_Subd subd)
	{
		return subd.maxDepth;
	}


	/*******************************************************************************
	 * VertexCount -- Retrieve the number of vertices
	 *
	 */
	static int ccs__VertexCount(cc_Mesh cage, int maxDepth)
	{
		return ccm_VertexCountAtDepth_Fast(cage, maxDepth);
	}

	public static int ccs_VertexCount(cc_Subd subd)
	{
		return ccs__VertexCount(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * Create -- Create a subd
	 *
	 */
	public static cc_Subd ccs_Create(int maxDepth, cc_Mesh cage)
	{
		int halfedgeCount = ccs_CumulativeHalfedgeCountAtDepth(cage,
																	maxDepth);
		int creaseCount = ccs_CumulativeCreaseCountAtDepth(cage,
																maxDepth);
		int vertexCount = ccs__VertexCount(cage, maxDepth);

		cc_Subd subd = new cc_Subd();
		subd.maxDepth = maxDepth;
		subd.halfedges = new cc_Halfedge_Quad[halfedgeCount];
		subd.creases = new cc_Crease[creaseCount];
		subd.vertexPoints = new float3[vertexCount];
		subd.cage = cage;

		return subd;
	}


	/*******************************************************************************
	 * Crease data accessors
	 *
	 * These accessors are hidden from the user because not all edges within
	 * the subd map to a crease. In particular: any edge create within a face
	 * does not have an associated crease. This is because such edges will never
	 * be sharp by construction.
	 *
	 */
	public static cc_Crease ccs__Crease(cc_Subd subd, int edgeID, int depth)
	{
		Debug.Assert(depth <= ccs_MaxDepth(subd) && depth > 0);
		int stride = ccs_CumulativeCreaseCountAtDepth(subd.cage,
																depth - 1);

		return subd.creases[stride + edgeID];
	}

	public static float ccs_CreaseSharpness_Fast(cc_Subd subd, int edgeID, int depth)
	{
		return ccs__Crease(subd, edgeID, depth).sharpness;
	}

	public static float ccs_CreaseSharpness(cc_Subd subd, int edgeID, int depth)
	{
		int creaseCount = ccm_CreaseCountAtDepth(subd.cage, depth);

		if (edgeID < creaseCount)
		{
			return ccs_CreaseSharpness_Fast(subd, edgeID, depth);
		}
		else
		{
			return 0.0f;
		}
	}

	public static int ccs_CreaseNextID_Fast(cc_Subd subd, int edgeID, int depth)
	{
		return ccs__Crease(subd, edgeID, depth).nextID;
	}

	public static int ccs_CreaseNextID(cc_Subd subd, int edgeID, int depth)
	{
		int creaseCount = ccm_CreaseCountAtDepth(subd.cage, depth);

		if (edgeID < creaseCount)
		{
			return ccs_CreaseNextID_Fast(subd, edgeID, depth);
		}
		else
		{
			return edgeID;
		}
	}

	public static int ccs_CreasePrevID_Fast(cc_Subd subd, int edgeID, int depth)
	{
		return ccs__Crease(subd, edgeID, depth).prevID;
	}

	public static int ccs_CreasePrevID(cc_Subd subd, int edgeID, int depth)
	{
		int creaseCount = ccm_CreaseCountAtDepth(subd.cage, depth);

		if (edgeID < creaseCount)
		{
			return ccs_CreasePrevID_Fast(subd, edgeID, depth);
		}
		else
		{
			return edgeID;
		}
	}


	/*******************************************************************************
	 * Halfedge data accessors
	 *
	 */
	public static cc_Halfedge_Quad ccs__Halfedge(cc_Subd subd, int halfedgeID, int depth)
	{
		Debug.Assert(depth <= ccs_MaxDepth(subd) && depth > 0);
		int stride = ccs_CumulativeHalfedgeCountAtDepth(subd.cage,
																  depth - 1);

		return subd.halfedges[stride + halfedgeID];
	}

	public static int ccs_HalfedgeVertexPointID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccs__Halfedge(subd, halfedgeID, depth).vertexID;
	}

	public static int ccs_HalfedgeTwinID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccs__Halfedge(subd, halfedgeID, depth).twinID;
	}

	public static int ccs_HalfedgeNextID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccm_HalfedgeNextID_Quad(halfedgeID);
	}

	public static int ccs_HalfedgePrevID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccm_HalfedgePrevID_Quad(halfedgeID);
	}

	public static int ccs_HalfedgeFaceID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccm_HalfedgeFaceID_Quad(halfedgeID);
	}

	public static int ccs_HalfedgeEdgeID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccs__Halfedge(subd, halfedgeID, depth).edgeID;
	}

	public static float ccs_HalfedgeSharpness(cc_Subd subd, int halfedgeID, int depth)
	{
		int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);

		return ccs_CreaseSharpness(subd, edgeID, depth);
	}

	public static float3 ccs_HalfedgeVertexPoint(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccs_VertexPoint(subd, ccs_HalfedgeVertexPointID(subd, halfedgeID, depth));
	}

	public static int ccs__HalfedgeVertexUvID(cc_Subd subd, int halfedgeID, int depth)
	{
		return ccs__Halfedge(subd, halfedgeID, depth).uvID;
	}

	public static float2 ccs_HalfedgeVertexUv(cc_Subd subd, int halfedgeID, int depth)
	{
		return cc__DecodeUv(ccs__Halfedge(subd, halfedgeID, depth).uvID);
	}


	/*******************************************************************************
	 * Vertex data accessors
	 *
	 */
	public static float3 ccs_VertexPoint(cc_Subd subd, int vertexID)
	{
		return subd.vertexPoints[vertexID];
	}


	/*******************************************************************************
	 * Vertex halfedge iteration
	 *
	 */
	public static int ccs_PrevVertexHalfedgeID(cc_Subd subd, int halfedgeID, int depth)
	{
		int prevID = ccs_HalfedgePrevID(subd, halfedgeID, depth);

		return ccs_HalfedgeTwinID(subd, prevID, depth);
	}

	public static int ccs_NextVertexHalfedgeID(cc_Subd subd, int halfedgeID, int depth)
	{
		int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);

		return ccs_HalfedgeNextID(subd, twinID, depth);
	}


	/*******************************************************************************
	 * Face to Halfedge Mapping
	 *
	 */
	public static int ccs_FaceToHalfedgeID(cc_Subd subd, int faceID, int depth)
	{
		return ccm_FaceToHalfedgeID_Quad(faceID);
	}


	/*******************************************************************************
	 * Edge to Halfedge Mapping
	 *
	 * This procedure returns one of the ID of one of the halfedge that constitutes
	 * the edge. This routine has O(depth) complexity.
	 *
	 */
	public static int ccs__EdgeToHalfedgeID_First(cc_Mesh cage, int edgeID)
	{
		int edgeCount = ccm_EdgeCount(cage);

		if /* [2E, 2E + H) */ (edgeID >= 2 * edgeCount)
		{
			int halfedgeID = edgeID - 2 * edgeCount;
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);

			return math.max(4 * halfedgeID + 1, 4 * nextID + 2);

		}
		else if /* */ ((edgeID & 1) == 1)
		{
			int halfedgeID = ccm_EdgeToHalfedgeID(cage, edgeID >> 1);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);

			return 4 * nextID + 3;

		}
		else /* */
		{
			int halfedgeID = ccm_EdgeToHalfedgeID(cage, edgeID >> 1);

			return 4 * halfedgeID + 0;
		}
	}

	public static int ccs_EdgeToHalfedgeID(cc_Subd subd, int edgeID, int depth)
	{
#if false // recursive version
		if (depth > 1) {
			int edgeCount = ccm_EdgeCountAtDepth_Fast(subd.cage, depth - 1);

			if /* [2E, 2E + H) */ (edgeID >= 2 * edgeCount) {
				int halfedgeID = edgeID - 2 * edgeCount;
				int nextID = ccm_NextFaceHalfedgeID_Quad(halfedgeID);

				return cc__Max(4 * halfedgeID + 1, 4 * nextID + 2);

			} else if /* [E, 2E) */ (edgeID >= edgeCount) {
				int halfedgeID = ccs_EdgeToHalfedgeID(subd,
														  edgeID >> 1,
														  depth - 1);
				int nextID = ccm_NextFaceHalfedgeID_Quad(halfedgeID);

				return 4 * nextID + 3;

			} else /* [0, E) */ {
				int halfedgeID = ccs_EdgeToHalfedgeID(subd, edgeID >> 1, depth - 1);

				return 4 * halfedgeID + 0;
			}
		} else {
			return ccs__EdgeToHalfedgeID_First(subd.cage, edgeID);
		}
#else // non-recursive version
		uint heap = 1u;
		int edgeHalfedgeID = 0;
		int heapDepth = depth;

		// build heap
		for (; heapDepth > 1; --heapDepth)
		{
			int edgeCount = ccm_EdgeCountAtDepth_Fast(subd.cage,
																heapDepth - 1);

			if /* [2E, 2E + H) */ (edgeID >= 2 * edgeCount)
			{
				int halfedgeID = edgeID - 2 * edgeCount;
				int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);

				edgeHalfedgeID = math.max(4 * halfedgeID + 1, 4 * nextID + 2);
				break;
			}
			else
			{
				heap = (heap << 1) | ((uint)edgeID & 1); // TODO: CHECK INT->UINT CONVERSION
				edgeID >>= 1;
			}
		}

		// initialize root cfg
		if (heapDepth == 1)
		{
			edgeHalfedgeID = ccs__EdgeToHalfedgeID_First(subd.cage, edgeID);
		}

		// read heap
		while (heap > 1u)
		{
			if ((heap & 1u) == 1u)
			{
				int nextID = ccm_HalfedgeNextID_Quad(edgeHalfedgeID);

				edgeHalfedgeID = 4 * nextID + 3;
			}
			else
			{
				edgeHalfedgeID = 4 * edgeHalfedgeID + 0;
			}

			heap >>= 1;
		}

		return edgeHalfedgeID;
#endif
	}


	/*******************************************************************************
	 * Vertex to Halfedge Mapping
	 *
	 * This procedure returns the ID of one of the halfedge that connects a
	 * given vertex. This routine has O(depth) complexity.
	 *
	 */
	public static int ccs__VertexToHalfedgeID_First(cc_Mesh cage, int vertexID)
	{
		int vertexCount = ccm_VertexCount(cage);
		int faceCount = ccm_FaceCount(cage);

		if /* [V + F, V + F + E) */ (vertexID >= vertexCount + faceCount)
		{
			int edgeID = vertexID - vertexCount - faceCount;

			return 4 * ccm_EdgeToHalfedgeID(cage, edgeID) + 1;

		}
		else if /* [V, V + F) */ (vertexID >= vertexCount)
		{
			int faceID = vertexID - vertexCount;

			return 4 * ccm_FaceToHalfedgeID(cage, faceID) + 2;

		}
		else /* [0, V) */
		{

			return 4 * ccm_VertexPointToHalfedgeID(cage, vertexID) + 0;
		}
	}

	public static int ccs_VertexPointToHalfedgeID(cc_Subd subd, int vertexID, int depth)
	{
#if false // recursive version
		if (depth > 1) {
			cc_Mesh cage = subd.cage;
			int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth - 1);
			int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth - 1);

			if /* [V + F, V + F + E) */ (vertexID >= vertexCount + faceCount) {
				int edgeID = vertexID - vertexCount - faceCount;

				return 4 * ccs_EdgeToHalfedgeID(subd, edgeID, depth - 1) + 1;

			} else if /* [V, V + F) */ (vertexID >= vertexCount) {
				int faceID = vertexID - vertexCount;

				return 4 * ccm_FaceToHalfedgeID_Quad(faceID) + 2;

			} else /* [0, V) */ {

				return 4 * ccs_VertexPointToHalfedgeID(subd, vertexID, depth - 1) + 0;
			}
		} else {

			return ccs__VertexToHalfedgeID_First(subd.cage, vertexID);
		}
#else // non-recursive version
		cc_Mesh cage = subd.cage;
		int heapDepth = depth;
		int stride = 0;
		int halfedgeID = 0; // TODO: CHECK IF CORRECT

		// build heap
		for (; heapDepth > 1; --heapDepth)
		{
			int vertexCount = ccm_VertexCountAtDepth_Fast(cage, heapDepth - 1);
			int faceCount = ccm_FaceCountAtDepth_Fast(cage, heapDepth - 1);

			if /* [V + F, V + F + E) */ (vertexID >= vertexCount + faceCount)
			{
				int edgeID = vertexID - faceCount - vertexCount;

				halfedgeID = 4 * ccs_EdgeToHalfedgeID(subd, edgeID, heapDepth - 1) + 1;
				break;
			}
			else if /* [V, V + F) */ (vertexID >= vertexCount)
			{
				int faceID = vertexID - vertexCount;

				halfedgeID = 4 * ccm_FaceToHalfedgeID_Quad(faceID) + 2;
				break;
			}
			else /* [0, V) */
			{
				stride += 2;
			}
		}

		// initialize root cfg
		if (heapDepth == 1)
		{
			halfedgeID = ccs__VertexToHalfedgeID_First(subd.cage, vertexID);
		}

		return halfedgeID << stride;
#endif
	}


	/*******************************************************************************
	 * CageFacePoints -- Applies Catmull Clark's face rule on the cage mesh
	 *
	 * The "Gather" routine iterates over each face of the mesh and compute the
	 * resulting face vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the face vertex.
	 *
	 */
	public static void ccs__CageFacePoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int faceID = 0; faceID < faceCount; ++faceID)
		{
			int halfedgeID = ccm_FaceToHalfedgeID(cage, faceID);
			float3 newFacePoint = ccm_HalfedgeVertexPoint(cage, halfedgeID);
			float faceVertexCount = 1.0f;

			for (int halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeID);
						 halfedgeIt != halfedgeID;
						 halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeIt))
			{
				float3 vertexPoint = ccm_HalfedgeVertexPoint(cage, halfedgeIt);

				newFacePoint += vertexPoint;

				++faceVertexCount;
			}

			newFacePoint /= faceVertexCount;

			subd.vertexPoints[newFacePointsStart + faceID] = newFacePoint;
		}
		//CC_BARRIER
	}

	public static void ccs__CageFacePoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			float3 vertexPoint = ccm_HalfedgeVertexPoint(cage, halfedgeID);
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			float faceVertexCount = 1.0f;

			for (int halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeID);
						 halfedgeIt != halfedgeID;
						 halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeIt))
			{
				++faceVertexCount;
			}

			//CC_ATOMIC
			subd.vertexPoints[newFacePointsStart + faceID] += vertexPoint / (float)faceVertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * CageEdgePoints -- Applies Catmull Clark's edge rule on the cage mesh
	 *
	 * The "Gather" routine iterates over each edge of the mesh and computes the
	 * resulting edge vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the edge vertex.
	 *
	 */
	public static void ccs__CageEdgePoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int edgeCount = ccm_EdgeCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int halfedgeID = ccm_EdgeToHalfedgeID(cage, edgeID);
			int twinID = ccm_HalfedgeTwinID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			float edgeWeight = twinID < 0 ? 0.0f : 1.0f;
			float3[] oldEdgePoints = new float3[]
			{
				ccm_HalfedgeVertexPoint(cage, halfedgeID),
				ccm_HalfedgeVertexPoint(cage, nextID)
			};
			float3[] newFacePointPair = new float3[]
			{
				subd.vertexPoints[newFacePointsStart + ccm_HalfedgeFaceID(cage, halfedgeID)],
				subd.vertexPoints[newFacePointsStart + ccm_HalfedgeFaceID(cage, math.max(0, twinID))]
			};
			float3 sharpEdgePoint = float3.zero;
			float3 smoothEdgePoint = float3.zero;
			float3 tmp1, tmp2;

			tmp1 = oldEdgePoints[0] + oldEdgePoints[1];
			tmp2 = newFacePointPair[0] + newFacePointPair[1];
			sharpEdgePoint = tmp1 * 0.5f;
			smoothEdgePoint = tmp1 + tmp2;
			smoothEdgePoint = smoothEdgePoint * 0.25f;
			subd.vertexPoints[newEdgePointsStart + edgeID] = math.lerp(sharpEdgePoint, smoothEdgePoint, edgeWeight);
		}

		//CC_BARRIER
	}

	public static void ccs__CageEdgePoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int faceCount = ccm_FaceCount(cage);
		int vertexCount = ccm_VertexCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int twinID = ccm_HalfedgeTwinID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3 tmp1, tmp2, tmp3, tmp4, atomicWeight;
			float weight = twinID >= 0 ? 0.5f : 1.0f;

			tmp1 = newFacePoint * 0.5f;
			tmp2 = ccm_HalfedgeVertexPoint(cage, halfedgeID) * weight;
			tmp3 = ccm_HalfedgeVertexPoint(cage, nextID) * weight;
			tmp4 = math.lerp(tmp2, tmp3, 0.5f);
			atomicWeight = math.lerp(tmp1, tmp4, weight);

			//CC_ATOMIC
			subd.vertexPoints[newEdgePointsStart + edgeID] += atomicWeight;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * CageEdgePoints -- Applies DeRole et al.'s edge rule on the cage mesh
	 *
	 * The "Gather" routine iterates over each edge of the mesh and computes the
	 * resulting edge vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the edge vertex.
	 *
	 */
	public static void ccs__CreasedCageEdgePoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int edgeCount = ccm_EdgeCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int halfedgeID = ccm_EdgeToHalfedgeID(cage, edgeID);
			int twinID = ccm_HalfedgeTwinID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			float sharp = ccm_CreaseSharpness(cage, edgeID);
			float edgeWeight = math.saturate(sharp);
			float3[] oldEdgePoints = new float3[]
			{
				ccm_HalfedgeVertexPoint(cage, halfedgeID),
				ccm_HalfedgeVertexPoint(cage,     nextID)
			};
			float3[] newAdjacentFacePoints = new float3[]
			{
				subd.vertexPoints[newFacePointsStart + ccm_HalfedgeFaceID(cage, halfedgeID)],
				subd.vertexPoints[newFacePointsStart + ccm_HalfedgeFaceID(cage, math.max(0, twinID))]
			};
			float3 sharpEdgePoint = float3.zero;
			float3 smoothEdgePoint = float3.zero;
			float3 tmp1, tmp2;

			tmp1 = oldEdgePoints[0] + oldEdgePoints[1];
			tmp2 = newAdjacentFacePoints[0] + newAdjacentFacePoints[1];
			sharpEdgePoint = tmp1 * 0.5f;
			smoothEdgePoint = tmp1 + tmp2;
			smoothEdgePoint = smoothEdgePoint * 0.25f;
			subd.vertexPoints[newEdgePointsStart + edgeID] = math.lerp(smoothEdgePoint, sharpEdgePoint, edgeWeight);
		}
		//CC_BARRIER
	}

	public static void ccs__CreasedCageEdgePoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int faceCount = ccm_FaceCount(cage);
		int vertexCount = ccm_VertexCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int twinID = ccm_HalfedgeTwinID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			float sharp = ccm_CreaseSharpness(cage, edgeID);
			float edgeWeight = math.saturate(sharp);
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3[] oldEdgePoints = new float3[]
			{
				ccm_HalfedgeVertexPoint(cage, halfedgeID),
				ccm_HalfedgeVertexPoint(cage, nextID)
			};
			float3 smoothPoint = float3.zero;
			float3 sharpPoint = float3.zero;
			float3 tmp, atomicWeight;

			// sharp point
			tmp = math.lerp(oldEdgePoints[0], oldEdgePoints[1], 0.5f);
			sharpPoint = tmp * (twinID < 0 ? 1.0f : 0.5f);

			// smooth point
			tmp = math.lerp(oldEdgePoints[0], newFacePoint, 0.5f);
			smoothPoint = tmp * 0.5f;

			// atomic weight
			atomicWeight = math.lerp(
				smoothPoint,
				sharpPoint,
				edgeWeight);

			//CC_ATOMIC
			subd.vertexPoints[newEdgePointsStart + edgeID] += atomicWeight;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * CageVertexPoints -- Applies Catmull Clark's vertex rule on the cage mesh
	 *
	 * The "Gather" routine iterates over each vertex of the mesh and computes the
	 * resulting smooth vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the smooth vertex.
	 *
	 */
	public static void ccs__CageVertexPoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int vertexID = 0; vertexID < vertexCount; ++vertexID)
		{
			int halfedgeID = ccm_VertexPointToHalfedgeID(cage, vertexID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			float3 newEdgePoint = subd.vertexPoints[newEdgePointsStart + edgeID];
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3 oldVertexPoint = cage.vertexPoints[vertexID];
			float3 smoothPoint = float3.zero;
			float valence = 1.0f;
			int iterator;
			float3 tmp1, tmp2;

			tmp1 = newFacePoint * -1.0f;
			tmp2 = newEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;

			for (iterator = ccm_PrevVertexHalfedgeID(cage, halfedgeID);
				 iterator >= 0 && iterator != halfedgeID;
				 iterator = ccm_PrevVertexHalfedgeID(cage, iterator))
			{
				int edgeID2 = ccm_HalfedgeEdgeID(cage, iterator);
				int faceID2 = ccm_HalfedgeFaceID(cage, iterator);
				int newEdgePointVertexID = vertexCount + faceCount + edgeID2;
				int newFacePointVertexID = vertexCount + faceID2;
				float3 newEdgePoint2 = ccs_VertexPoint(subd, newEdgePointVertexID);
				float3 newFacePoint2 = ccs_VertexPoint(subd, newFacePointVertexID);

				tmp1 = newFacePoint2 * -1.0f;
				tmp2 = newEdgePoint2 * 4.0f;
				smoothPoint = smoothPoint + tmp1;
				smoothPoint = smoothPoint + tmp2;
				++valence;
			}

			tmp1 = smoothPoint * 1.0f / (valence * valence);
			tmp2 = oldVertexPoint * 1.0f - 3.0f / valence;
			smoothPoint = tmp1 + tmp2;
			subd.vertexPoints[vertexID] = math.lerp(
				oldVertexPoint,
				smoothPoint,
				iterator != halfedgeID ? 0.0f : 1.0f);
		}
		//CC_BARRIER
	}

	public static void ccs__CageVertexPoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int faceCount = ccm_FaceCount(cage);
		int vertexCount = ccm_VertexCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccm_HalfedgeVertexID(cage, halfedgeID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			int valence = 1;
			int forwardIterator, backwardIterator;

			for (forwardIterator = ccm_PrevVertexHalfedgeID(cage, halfedgeID);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccm_PrevVertexHalfedgeID(cage, forwardIterator))
			{
				++valence;
			}

			for (backwardIterator = ccm_NextVertexHalfedgeID(cage, halfedgeID);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccm_NextVertexHalfedgeID(cage, backwardIterator))
			{
				++valence;
			}

			float w = 1.0f / (float)valence;
			float3 v = cage.vertexPoints[vertexID];
			float3 f = subd.vertexPoints[newFacePointsStart + faceID];
			float3 e = subd.vertexPoints[newEdgePointsStart + edgeID];
			float s = forwardIterator < 0 ? 0.0f : 1.0f;
			//CC_ATOMIC
			subd.vertexPoints[vertexID] +=
				w * (v + w * s * (4.0f * e - f - 3.0f * v));

		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * CreasedCageVertexPoints -- Applies DeRose et al.'s vertex rule on cage mesh
	 *
	 * The "Gather" routine iterates over each vertex of the mesh and computes the
	 * resulting smooth vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the smooth vertex.
	 *
	 */
	public static void ccs__CreasedCageVertexPoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int vertexID = 0; vertexID < vertexCount; ++vertexID)
		{
			int halfedgeID = ccm_VertexPointToHalfedgeID(cage, vertexID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int prevID = ccm_HalfedgePrevID(cage, halfedgeID);
			int prevEdgeID = ccm_HalfedgeEdgeID(cage, prevID);
			int prevFaceID = ccm_HalfedgeFaceID(cage, prevID);
			float thisS = ccm_HalfedgeSharpness(cage, halfedgeID);
			float prevS = ccm_HalfedgeSharpness(cage, prevID);
			float creaseWeight = math.sign(thisS);
			float prevCreaseWeight = math.sign(prevS);
			float3 newEdgePoint = subd.vertexPoints[newEdgePointsStart + edgeID];
			float3 newPrevEdgePoint = subd.vertexPoints[newEdgePointsStart + prevEdgeID];
			float3 newPrevFacePoint = subd.vertexPoints[newFacePointsStart + prevFaceID];
			float3 oldPoint = ccm_VertexPoint(cage, vertexID);
			float3 smoothPoint = float3.zero;
			float3 creasePoint = float3.zero;
			float avgS = prevS;
			float creaseCount = prevCreaseWeight;
			float valence = 1.0f;
			int forwardIterator;
			float3 tmp1, tmp2;

			// smooth contrib
			tmp1 = newPrevFacePoint * -1.0f;
			tmp2 = newPrevEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;

			// crease contrib
			tmp1 = newPrevEdgePoint * prevCreaseWeight;
			creasePoint = creasePoint + tmp1;

			for (forwardIterator = ccm_HalfedgeTwinID(cage, prevID);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccm_HalfedgeTwinID(cage, forwardIterator))
			{
				int prevID2 = ccm_HalfedgePrevID(cage, forwardIterator);
				int prevEdgeID2 = ccm_HalfedgeEdgeID(cage, prevID2);
				int prevFaceID2 = ccm_HalfedgeFaceID(cage, prevID2);
				float3 newPrevEdgePoint2 = subd.vertexPoints[newEdgePointsStart + prevEdgeID2];
				float3 newPrevFacePoint2 = subd.vertexPoints[newFacePointsStart + prevFaceID2];
				float prevS2 = ccm_HalfedgeSharpness(cage, prevID2);
				float prevCreaseWeight2 = math.sign(prevS2);

				// smooth contrib
				tmp1 = newPrevFacePoint2 * -1.0f;
				tmp2 = newPrevEdgePoint2 * 4.0f;
				smoothPoint = smoothPoint + tmp1;
				smoothPoint = smoothPoint + tmp2;
				++valence;

				// crease contrib
				tmp1 = newPrevEdgePoint2 * prevCreaseWeight2;
				creasePoint = creasePoint + tmp1;
				avgS += prevS2;
				creaseCount += prevCreaseWeight2;

				// next vertex halfedge
				forwardIterator = prevID2;
			}

			// boundary corrections
			if (forwardIterator < 0)
			{
				tmp1 = newEdgePoint * creaseWeight;
				creasePoint = creasePoint + tmp1;
				creaseCount += creaseWeight;
				++valence;
			}

			// average sharpness
			avgS /= valence;

			// smooth point
			tmp1 = smoothPoint * 1.0f / (valence * valence);
			tmp2 = oldPoint * (1.0f - 3.0f / valence);
			smoothPoint = tmp1 + tmp2;

			// crease point
			tmp1 = creasePoint * 0.25f;
			tmp2 = oldPoint * 0.5f;
			creasePoint = tmp1 + tmp2;

			// proper vertex rule selection (TODO: make branchless)
			if (creaseCount <= 1.0f)
			{
				subd.vertexPoints[vertexID] = smoothPoint;
			}
			else if (creaseCount >= 3.0f || valence == 2.0f)
			{
				subd.vertexPoints[vertexID] = oldPoint;
			}
			else
			{
				subd.vertexPoints[vertexID] = math.lerp(
					oldPoint,
					creasePoint,
					math.saturate(avgS));
			}
		}
		//CC_BARRIER
	}

	public static void ccs__CreasedCageVertexPoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int faceCount = ccm_FaceCount(cage);
		int vertexCount = ccm_VertexCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccm_HalfedgeVertexID(cage, halfedgeID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			int prevID = ccm_HalfedgePrevID(cage, halfedgeID);
			int prevEdgeID = ccm_HalfedgeEdgeID(cage, prevID);
			float thisS = ccm_HalfedgeSharpness(cage, halfedgeID);
			float prevS = ccm_HalfedgeSharpness(cage, prevID);
			float creaseWeight = math.sign(thisS);
			float prevCreaseWeight = math.sign(prevS);
			float3 newPrevEdgePoint = subd.vertexPoints[newEdgePointsStart + prevEdgeID];
			float3 newEdgePoint = subd.vertexPoints[newEdgePointsStart + edgeID];
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3 oldPoint = cage.vertexPoints[vertexID];
			float3 cornerPoint = float3.zero;
			float3 smoothPoint = float3.zero;
			float3 creasePoint = float3.zero;
			float3 atomicWeight = float3.zero;
			float avgS = prevS;
			float creaseCount = prevCreaseWeight;
			float valence = 1.0f;
			int forwardIterator, backwardIterator;
			float3 tmp1, tmp2;

			for (forwardIterator = ccm_HalfedgeTwinID(cage, prevID);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccm_HalfedgeTwinID(cage, forwardIterator))
			{
				int prevID2 = ccm_HalfedgePrevID(cage, forwardIterator);
				float prevS2 = ccm_HalfedgeSharpness(cage, prevID2);
				float prevCreaseWeight2 = math.sign(prevS2);

				// valence computation
				++valence;

				// crease computation
				avgS += prevS2;
				creaseCount += prevCreaseWeight2;

				// next vertex halfedge
				forwardIterator = prevID2;
			}

			for (backwardIterator = ccm_HalfedgeTwinID(cage, halfedgeID);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccm_HalfedgeTwinID(cage, backwardIterator))
			{
				int nextID = ccm_HalfedgeNextID(cage, backwardIterator);
				float nextS = ccm_HalfedgeSharpness(cage, nextID);
				float nextCreaseWeight = math.sign(nextS);

				// valence computation
				++valence;

				// crease computation
				avgS += nextS;
				creaseCount += nextCreaseWeight;

				// next vertex halfedge
				backwardIterator = nextID;
			}

			// corner point
			cornerPoint = oldPoint * 1.0f / valence;

			// crease computation: V / 4
			tmp1 = oldPoint * 0.25f * creaseWeight;
			tmp2 = newEdgePoint * 0.25f * creaseWeight;
			creasePoint = tmp1 + tmp2;

			// smooth computation: (4E - F + (n - 3) V) / N
			tmp1 = newFacePoint * -1.0f;
			tmp2 = newEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;
			tmp1 = oldPoint * valence - 3.0f;
			smoothPoint = smoothPoint + tmp1;
			smoothPoint = smoothPoint * 1.0f / (valence * valence);

			// boundary corrections
			if (forwardIterator < 0)
			{
				creaseCount += creaseWeight;
				++valence;

				tmp1 = oldPoint * 0.25f * prevCreaseWeight;
				tmp2 = newPrevEdgePoint * 0.25f * prevCreaseWeight;
				tmp1 = tmp1 + tmp2;
				creasePoint = creasePoint + tmp1;
			}

			// average sharpness
			avgS /= valence;

			// atomicWeight (TODO: make branchless ?)
			if (creaseCount <= 1.0f)
			{
				atomicWeight = smoothPoint;
			}
			else if (creaseCount >= 3.0f || valence == 2.0f)
			{
				atomicWeight = cornerPoint;
			}
			else
			{
				atomicWeight = math.lerp(
					cornerPoint,
					creasePoint,
					math.saturate(avgS));
			}

			//CC_ATOMIC
			subd.vertexPoints[vertexID] += atomicWeight;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * FacePoints -- Applies Catmull Clark's face rule on the subd
	 *
	 * The "Gather" routine iterates over each face of the mesh and compute the
	 * resulting face vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the face vertex.
	 *
	 */
	public static void ccs__FacePoints_Gather(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int faceID = 0; faceID < faceCount; ++faceID)
		{
			int faceVertexID = vertexCount + faceID;
			subd.vertexPoints[faceVertexID] = float3.zero;

			for (int faceVertexID2 = 0; faceVertexID2 < 4; ++faceVertexID2)
			{
				int halfedgeID = 4 * faceID + faceVertexID2;
				int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);

				subd.vertexPoints[faceVertexID] += subd.vertexPoints[vertexID] / 4.0f;
			}
		}
		//CC_BARRIER
	}

	public static void ccs__FacePoints_Scatter(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int faceID = ccm_HalfedgeFaceID_Quad(halfedgeID);
			int faceVertexID = vertexCount + faceID;

			//CC_ATOMIC
			subd.vertexPoints[faceVertexID] += subd.vertexPoints[vertexID] / 4.0f;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * EdgePoints -- Applies Catmull Clark's edge rule on the subd
	 *
	 * The "Gather" routine iterates over each edge of the mesh and compute the
	 * resulting edge vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the edge vertex.
	 *
	 */
	public static void ccs__EdgePoints_Gather(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int edgeCount = ccm_EdgeCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int halfedgeID = ccs_EdgeToHalfedgeID(subd, edgeID, depth);
			int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
			int nextID = ccs_HalfedgeNextID(subd, halfedgeID, depth);
			float edgeWeight = twinID < 0 ? 0.0f : 1.0f;
			float3[] oldEdgePoints = new float3[] {
				ccs_HalfedgeVertexPoint(subd, halfedgeID, depth),
				ccs_HalfedgeVertexPoint(subd,     nextID, depth)
			};
			float3[] newAdjacentFacePoints = new float3[] {
				subd.vertexPoints[newFacePointsStart + ccs_HalfedgeFaceID(subd,          halfedgeID, depth)],
				subd.vertexPoints[newFacePointsStart + ccs_HalfedgeFaceID(subd, math.max(0, twinID), depth)]
			};
			float3 sharpEdgePoint = float3.zero;
			float3 smoothEdgePoint = float3.zero;
			float3 tmp1, tmp2;

			tmp1 = oldEdgePoints[0] + oldEdgePoints[1];
			tmp2 = newAdjacentFacePoints[0] + newAdjacentFacePoints[1];
			sharpEdgePoint = tmp1 * 0.5f;
			smoothEdgePoint = tmp1 + tmp2;
			smoothEdgePoint = smoothEdgePoint * 0.25f;
			subd.vertexPoints[newEdgePointsStart + edgeID] = math.lerp(
				sharpEdgePoint,
				smoothEdgePoint,
				edgeWeight);
		}
		//CC_BARRIER
	}

	public static void ccs__EdgePoints_Scatter(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int faceID = ccm_HalfedgeFaceID_Quad(halfedgeID);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int faceVertexID = vertexCount + faceID;
			int edgeVertexID = vertexCount + faceCount + edgeID;
			int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
			int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);
			float3 faceVertex = subd.vertexPoints[faceVertexID];
			float3[] edgeVertices = new float3[] {
				subd.vertexPoints[ccs_HalfedgeVertexPointID(subd, halfedgeID, depth)],
				subd.vertexPoints[ccs_HalfedgeVertexPointID(subd, nextID, depth)]
			};
			float3 tmp1, tmp2, tmp3, tmp4, edgeVertexWeight;
			float weight = twinID >= 0 ? 0.5f : 1.0f;

			tmp1 = faceVertex * 0.5f;
			tmp2 = edgeVertices[0] * weight;
			tmp3 = edgeVertices[1] * weight;
			tmp4 = math.lerp(tmp2, tmp3, 0.5f);
			edgeVertexWeight = math.lerp(tmp1, tmp4, weight);

			//CC_ATOMIC
			subd.vertexPoints[edgeVertexID] += edgeVertexWeight;
		}
		//CC_BARRIER
	}

	/*******************************************************************************
	 * CreasedEdgePoints -- Applies DeRose et al's edge rule on the subd
	 *
	 * The "Gather" routine iterates over each edge of the mesh and compute the
	 * resulting edge vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the edge vertex.
	 *
	 */
	public static void ccs__CreasedEdgePoints_Gather(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int edgeCount = ccm_EdgeCountAtDepth_Fast(cage, depth);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int halfedgeID = ccs_EdgeToHalfedgeID(subd, edgeID, depth);
			int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
			int nextID = ccs_HalfedgeNextID(subd, halfedgeID, depth);
			float sharp = ccs_CreaseSharpness(subd, edgeID, depth);
			float edgeWeight = math.saturate(sharp);
			float3[] oldEdgePoints = new float3[] {
				ccs_HalfedgeVertexPoint(subd, halfedgeID, depth),
				ccs_HalfedgeVertexPoint(subd,     nextID, depth)
			};
			float3[] newAdjacentFacePoints = new float3[] {
				subd.vertexPoints[newFacePointsStart + ccs_HalfedgeFaceID(subd,          halfedgeID, depth)],
				subd.vertexPoints[newFacePointsStart + ccs_HalfedgeFaceID(subd, math.max(0, twinID), depth)]
			};
			float3 sharpEdgePoint = float3.zero;
			float3 smoothEdgePoint = float3.zero;
			float3 tmp1, tmp2;

			tmp1 = oldEdgePoints[0] + oldEdgePoints[1];
			tmp2 = newAdjacentFacePoints[0] + newAdjacentFacePoints[1];
			sharpEdgePoint = tmp1 * 0.5f;
			smoothEdgePoint = tmp1 + tmp2;
			smoothEdgePoint = smoothEdgePoint * 0.25f;
			subd.vertexPoints[newEdgePointsStart + edgeID] = math.lerp(
				smoothEdgePoint,
				sharpEdgePoint,
				edgeWeight);
		}
		//CC_BARRIER
	}

	public static void ccs__CreasedEdgePoints_Scatter(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int faceID = ccs_HalfedgeFaceID(subd, halfedgeID, depth);
			int nextID = ccs_HalfedgeNextID(subd, halfedgeID, depth);
			float sharp = ccs_CreaseSharpness(subd, edgeID, depth);
			float edgeWeight = math.saturate(sharp);
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3[] oldEdgePoints = new float3[] {
				ccs_HalfedgeVertexPoint(subd, halfedgeID, depth),
				ccs_HalfedgeVertexPoint(subd,     nextID, depth)
			};
			float3 smoothPoint = float3.zero;
			float3 sharpPoint = float3.zero;
			float3 tmp, atomicWeight;

			// sharp point
			tmp = math.lerp(oldEdgePoints[0], oldEdgePoints[1], 0.5f);
			sharpPoint = tmp * (twinID < 0 ? 1.0f : 0.5f);

			// smooth point
			tmp = math.lerp(oldEdgePoints[0], newFacePoint, 0.5f);
			smoothPoint = tmp * 0.5f;

			// atomic weight
			atomicWeight = math.lerp(
				smoothPoint,
				sharpPoint,
				edgeWeight);

			//CC_ATOMIC
			subd.vertexPoints[newEdgePointsStart + edgeID] += atomicWeight;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * VertexPoints -- Applies Catmull Clark's vertex rule on the subd
	 *
	 * The "Gather" routine iterates over each vertex of the mesh and computes the
	 * resulting smooth vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the smooth vertex.
	 *
	 */
	public static void ccs__VertexPoints_Gather(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int vertexID = 0; vertexID < vertexCount; ++vertexID)
		{
			int halfedgeID = ccs_VertexPointToHalfedgeID(subd, vertexID, depth);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int faceID = ccs_HalfedgeFaceID(subd, halfedgeID, depth);
			int newEdgePointVertexID = vertexCount + faceCount + edgeID;
			int newFacePointVertexID = vertexCount + faceID;
			float3 newEdgePoint = ccs_VertexPoint(subd, newEdgePointVertexID);
			float3 newFacePoint = ccs_VertexPoint(subd, newFacePointVertexID);
			float3 oldVertexPoint = ccs_VertexPoint(subd, vertexID);
			float3 smoothPoint = float3.zero;
			float valence = 1.0f;
			int iterator;
			float3 tmp1, tmp2;

			tmp1 = newFacePoint * -1.0f;
			tmp2 = newEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;

			for (iterator = ccs_PrevVertexHalfedgeID(subd, halfedgeID, depth);
				 iterator >= 0 && iterator != halfedgeID;
				 iterator = ccs_PrevVertexHalfedgeID(subd, iterator, depth))
			{
				int edgeID2 = ccs_HalfedgeEdgeID(subd, iterator, depth);
				int faceID2 = ccs_HalfedgeFaceID(subd, iterator, depth);
				int newEdgePointVertexID2 = vertexCount + faceCount + edgeID2;
				int newFacePointVertexID2 = vertexCount + faceID2;
				float3 newEdgePoint2 = ccs_VertexPoint(subd, newEdgePointVertexID2);
				float3 newFacePoint2 = ccs_VertexPoint(subd, newFacePointVertexID2);

				tmp1 = newFacePoint2 * -1.0f;
				tmp2 = newEdgePoint2 * 4.0f;
				smoothPoint = smoothPoint + tmp1;
				smoothPoint = smoothPoint + tmp2;
				++valence;
			}

			tmp1 = smoothPoint * 1.0f / (valence * valence);
			tmp2 = oldVertexPoint * 1.0f - 3.0f / valence;
			smoothPoint = tmp1 + tmp2;
			subd.vertexPoints[vertexID] = math.lerp(
				oldVertexPoint,
				smoothPoint,
				iterator != halfedgeID ? 0.0f : 1.0f);
		}
		//CC_BARRIER
	}

	public static void ccs__VertexPoints_Scatter(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int valence = 1;
			int forwardIterator, backwardIterator;

			for (forwardIterator = ccs_PrevVertexHalfedgeID(subd, halfedgeID, depth);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccs_PrevVertexHalfedgeID(subd, forwardIterator, depth))
			{
				++valence;
			}

			for (backwardIterator = ccs_NextVertexHalfedgeID(subd, halfedgeID, depth);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccs_NextVertexHalfedgeID(subd, backwardIterator, depth))
			{
				++valence;
			}

			if (forwardIterator >= 0 && valence > 2)
			{
				float w = 1.0f / (float)valence;

				//CC_ATOMIC
				subd.vertexPoints[vertexID] *= math.pow(math.max(0.0f, 1.0f - 3.0f * w), w);
			}
		}
		//CC_BARRIER

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int faceID = ccm_HalfedgeFaceID_Quad(halfedgeID);
			int faceVertexID = vertexCount + faceID;
			int edgeVertexID = vertexCount + faceCount + edgeID;
			float3 faceVertex = subd.vertexPoints[faceVertexID];
			float3 edgeVertex = subd.vertexPoints[edgeVertexID];
			int valence = 1;
			int forwardIterator, backwardIterator;

			for (forwardIterator = ccs_PrevVertexHalfedgeID(subd, halfedgeID, depth);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccs_PrevVertexHalfedgeID(subd, forwardIterator, depth))
			{
				++valence;
			}

			for (backwardIterator = ccs_NextVertexHalfedgeID(subd, halfedgeID, depth);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccs_NextVertexHalfedgeID(subd, backwardIterator, depth))
			{
				++valence;
			}

			if (forwardIterator >= 0 && valence > 2)
			{
				float w = 1.0f / (float)(valence * valence);
				float3 e = edgeVertex;
				float3 f = faceVertex;

				//CC_ATOMIC
				subd.vertexPoints[vertexID] += w * (4.0f * e - f);
			}
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * CreasedVertexPoints -- Applies DeRose et al.'s vertex rule on the subd
	 *
	 * The "Gather" routine iterates over each vertex of the mesh and computes the
	 * resulting smooth vertex.
	 *
	 * The "Scatter" routine iterates over each halfedge of the mesh and atomically
	 * adds its contribution to the computation of the smooth vertex.
	 *
	 */
	public static void ccs__CreasedVertexPoints_Gather(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int vertexID = 0; vertexID < vertexCount; ++vertexID)
		{
			int halfedgeID = ccs_VertexPointToHalfedgeID(subd, vertexID, depth);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int prevID = ccs_HalfedgePrevID(subd, halfedgeID, depth);
			int prevEdgeID = ccs_HalfedgeEdgeID(subd, prevID, depth);
			int prevFaceID = ccs_HalfedgeFaceID(subd, prevID, depth);
			float thisS = ccs_HalfedgeSharpness(subd, halfedgeID, depth);
			float prevS = ccs_HalfedgeSharpness(subd, prevID, depth);
			float creaseWeight = math.sign(thisS);
			float prevCreaseWeight = math.sign(prevS);
			float3 newEdgePoint = subd.vertexPoints[newEdgePointsStart + edgeID];
			float3 newPrevEdgePoint = subd.vertexPoints[newEdgePointsStart + prevEdgeID];
			float3 newPrevFacePoint = subd.vertexPoints[newFacePointsStart + prevFaceID];
			float3 oldPoint = ccs_VertexPoint(subd, vertexID);
			float3 smoothPoint = float3.zero;
			float3 creasePoint = float3.zero;
			float avgS = prevS;
			float creaseCount = prevCreaseWeight;
			float valence = 1.0f;
			int forwardIterator, backwardIterator;
			float3 tmp1, tmp2;

			// smooth contrib
			tmp1 = newPrevFacePoint * -1.0f;
			tmp2 = newPrevEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;

			// crease contrib
			tmp1 = newPrevEdgePoint * prevCreaseWeight;
			creasePoint = creasePoint + tmp1;

			for (forwardIterator = ccs_HalfedgeTwinID(subd, prevID, depth);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccs_HalfedgeTwinID(subd, forwardIterator, depth))
			{
				int prevID2 = ccs_HalfedgePrevID(subd, forwardIterator, depth);
				int prevEdgeID2 = ccs_HalfedgeEdgeID(subd, prevID2, depth);
				int prevFaceID2 = ccs_HalfedgeFaceID(subd, prevID2, depth);
				float3 newPrevEdgePoint2 = subd.vertexPoints[newEdgePointsStart + prevEdgeID2];
				float3 newPrevFacePoint2 = subd.vertexPoints[newFacePointsStart + prevFaceID2];
				float prevS2 = ccs_HalfedgeSharpness(subd, prevID2, depth);
				float prevCreaseWeight2 = math.sign(prevS2);

				// smooth contrib
				tmp1 = newPrevFacePoint2 * -1.0f;
				tmp2 = newPrevEdgePoint2 * 4.0f;
				smoothPoint = smoothPoint + tmp1;
				smoothPoint = smoothPoint + tmp2;
				++valence;

				// crease contrib
				tmp1 = newPrevEdgePoint2 * prevCreaseWeight2;
				creasePoint = creasePoint + tmp1;
				avgS += prevS2;
				creaseCount += prevCreaseWeight2;

				// next vertex halfedge
				forwardIterator = prevID2;
			}

			for (backwardIterator = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccs_HalfedgeTwinID(subd, backwardIterator, depth))
			{
				int nextID = ccs_HalfedgeNextID(subd, backwardIterator, depth);
				int nextEdgeID = ccs_HalfedgeEdgeID(subd, nextID, depth);
				int nextFaceID = ccs_HalfedgeFaceID(subd, nextID, depth);
				float3 newNextEdgePoint = subd.vertexPoints[newEdgePointsStart + nextEdgeID];
				float3 newNextFacePoint = subd.vertexPoints[newFacePointsStart + nextFaceID];
				float nextS = ccs_HalfedgeSharpness(subd, nextID, depth);
				float nextCreaseWeight = math.sign(nextS);

				// smooth contrib
				tmp1 = newNextFacePoint * -1.0f;
				tmp2 = newNextEdgePoint * 4.0f;
				smoothPoint = smoothPoint + tmp1;
				smoothPoint = smoothPoint + tmp2;
				++valence;

				// crease contrib
				tmp1 = newNextEdgePoint * nextCreaseWeight;
				creasePoint = creasePoint + tmp1;
				avgS += nextS;
				creaseCount += nextCreaseWeight;

				// next vertex halfedge
				backwardIterator = nextID;
			}

			// boundary corrections
			if (forwardIterator < 0)
			{
				tmp1 = newEdgePoint * creaseWeight;
				creasePoint = creasePoint + tmp1;
				creaseCount += creaseWeight;
				++valence;
			}

			// average sharpness
			avgS /= valence;

			// smooth point
			tmp1 = smoothPoint * (1.0f / (valence * valence));
			tmp2 = oldPoint * (1.0f - 3.0f / valence);
			smoothPoint = tmp1 + tmp2;

			// crease point
			tmp1 = creasePoint * 0.5f / creaseCount;
			tmp2 = oldPoint * 0.5f;
			creasePoint = tmp1 + tmp2;

			// proper vertex rule selection (TODO: make branchless)
			if (creaseCount <= 1.0f)
			{
				subd.vertexPoints[vertexID] = smoothPoint;
			}
			else if (creaseCount >= 3.0f || valence == 2.0f)
			{
				subd.vertexPoints[vertexID] = oldPoint;
			}
			else
			{
				subd.vertexPoints[vertexID] = math.lerp(
					oldPoint,
					creasePoint,
					math.saturate(avgS));
			}
		}
		//CC_BARRIER
	}

	public static void ccs__CreasedVertexPoints_Scatter(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int newFacePointsStart = vertexCount; // TODO: CHECK IF CHANGE IS CORRECT
		int newEdgePointsStart = vertexCount + faceCount; // TODO: CHECK IF CHANGE IS CORRECT

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int prevID = ccs_HalfedgePrevID(subd, halfedgeID, depth);
			float thisS = ccs_HalfedgeSharpness(subd, halfedgeID, depth);
			float prevS = ccs_HalfedgeSharpness(subd, prevID, depth);
			float creaseWeight = math.sign(thisS);
			float prevCreaseWeight = math.sign(prevS);
			float avgS = prevS;
			float creaseCount = prevCreaseWeight;
			float valence = 1.0f, valenceRcp;
			int forwardIterator, backwardIterator;
			float atomicWeight;

			for (forwardIterator = ccs_HalfedgeTwinID(subd, prevID, depth);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccs_HalfedgeTwinID(subd, forwardIterator, depth))
			{
				int prevID2 = ccs_HalfedgePrevID(subd, forwardIterator, depth);
				float prevS2 = ccs_HalfedgeSharpness(subd, prevID2, depth);
				float prevCreaseWeight2 = math.sign(prevS2);

				// valence computation
				++valence;

				// crease computation
				avgS += prevS2;
				creaseCount += prevCreaseWeight2;

				// next vertex halfedge
				forwardIterator = prevID2;
			}

			for (backwardIterator = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccs_HalfedgeTwinID(subd, backwardIterator, depth))
			{
				int nextID = ccs_HalfedgeNextID(subd, backwardIterator, depth);
				float nextS = ccs_HalfedgeSharpness(subd, nextID, depth);
				float nextCreaseWeight = math.sign(nextS);

				// valence computation
				++valence;

				// crease computation
				avgS += nextS;
				creaseCount += nextCreaseWeight;

				// next vertex halfedge
				backwardIterator = nextID;
			}

			// boundary corrections
			if (forwardIterator < 0)
			{
				creaseCount += creaseWeight;
				++valence;
			}

			// average sharpness
			valenceRcp = 1.0f / valence;
			avgS *= valenceRcp;

			// atomicWeight (TODO: make branchless ?)
			if (creaseCount <= 1.0f)
			{
				float w = valenceRcp;

				atomicWeight = math.pow(math.max(0.0f, 1.0f - 3.0f * w), w);
			}
			else if (creaseCount >= 3.0f || valence == 2.0f)
			{
				atomicWeight = 1.0f;
			}
			else
			{
				float w = math.saturate(avgS);
				float tmp = math.sqrt(1.0f - 0.5f * w);

				// TODO: cleanup !!!
				if (forwardIterator < 0 && creaseWeight > 0.0f && prevCreaseWeight > 0.0f)
					atomicWeight = 1.0f - 0.5f * w;
				else if ((forwardIterator < 0 && prevCreaseWeight > 0.0f) || creaseWeight > 0.0f)
					atomicWeight = tmp;
				else
					atomicWeight = 1.0f;
			}

			//CC_ATOMIC
			subd.vertexPoints[vertexID] *= atomicWeight;
		}
		//CC_BARRIER


		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int faceID = ccs_HalfedgeFaceID(subd, halfedgeID, depth);
			int prevID = ccs_HalfedgePrevID(subd, halfedgeID, depth);
			int prevEdgeID = ccs_HalfedgeEdgeID(subd, prevID, depth);
			float thisS = ccs_HalfedgeSharpness(subd, halfedgeID, depth);
			float prevS = ccs_HalfedgeSharpness(subd, prevID, depth);
			float creaseWeight = math.sign(thisS);
			float prevCreaseWeight = math.sign(prevS);
			float3 newPrevEdgePoint = subd.vertexPoints[newEdgePointsStart + prevEdgeID];
			float3 newEdgePoint = subd.vertexPoints[newEdgePointsStart + edgeID];
			float3 newFacePoint = subd.vertexPoints[newFacePointsStart + faceID];
			float3 smoothPoint = float3.zero;
			float3 creasePoint = float3.zero;
			float3 atomicWeight = float3.zero;
			float avgS = prevS;
			float creaseCount = prevCreaseWeight;
			float valence = 1.0f;
			int forwardIterator, backwardIterator;
			float3 tmp1, tmp2;

			for (forwardIterator = ccs_HalfedgeTwinID(subd, prevID, depth);
				 forwardIterator >= 0 && forwardIterator != halfedgeID;
				 forwardIterator = ccs_HalfedgeTwinID(subd, forwardIterator, depth))
			{
				int prevID2 = ccs_HalfedgePrevID(subd, forwardIterator, depth);
				float prevS2 = ccs_HalfedgeSharpness(subd, prevID2, depth);
				float prevCreaseWeight2 = math.sign(prevS2);

				// valence computation
				++valence;

				// crease computation
				avgS += prevS2;
				creaseCount += prevCreaseWeight2;

				// next vertex halfedge
				forwardIterator = prevID2;
			}

			for (backwardIterator = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
				 forwardIterator < 0 && backwardIterator >= 0 && backwardIterator != halfedgeID;
				 backwardIterator = ccs_HalfedgeTwinID(subd, backwardIterator, depth))
			{
				int nextID = ccs_HalfedgeNextID(subd, backwardIterator, depth);
				float nextS = ccs_HalfedgeSharpness(subd, nextID, depth);
				float nextCreaseWeight = math.sign(nextS);

				// valence computation
				++valence;

				// crease computation
				avgS += nextS;
				creaseCount += nextCreaseWeight;

				// next vertex halfedge
				backwardIterator = nextID;
			}

			// crease computation: V / 4
			creasePoint = newEdgePoint * 0.25f * creaseWeight;

			// smooth computation: (4E - F + (n - 3) V) / N
			tmp1 = newFacePoint * -1.0f;
			tmp2 = newEdgePoint * 4.0f;
			smoothPoint = tmp1 + tmp2;
			smoothPoint = smoothPoint * 1.0f / (valence * valence);

			// boundary corrections
			if (forwardIterator < 0)
			{
				creaseCount += creaseWeight;
				++valence;

				tmp1 = newPrevEdgePoint * 0.25f * prevCreaseWeight;
				creasePoint = creasePoint + tmp1;
			}

			// average sharpness
			avgS /= valence;

			// atomicWeight (TODO: make branchless ?)
			if (creaseCount <= 1.0f)
			{
				atomicWeight = smoothPoint;
			}
			else if (creaseCount >= 3.0f || valence == 2.0f)
			{
				// weight is already 0;
			}
			else
			{
				atomicWeight = creasePoint * math.saturate(avgS);
			}

			//CC_ATOMIC
			subd.vertexPoints[vertexID] += atomicWeight;
		}
		//CC_BARRIER
	}



	/*******************************************************************************
	 * RefineVertexPoints -- Computes the result of Catmull Clark subdivision.
	 *
	 */
	static void ccs__ClearVertexPoints(cc_Subd subd)
	{
		for (int i = 0; i < subd.vertexPoints.Length; i++)
			subd.vertexPoints[i] = float3.zero;
	}

	public static void ccs_RefineVertexPoints_Scatter(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;

		ccs__ClearVertexPoints(subd);

		if (true/*ccm_IsCreased(cage)*/)
		{
			ccs__CageFacePoints_Scatter(subd);
			ccs__CreasedCageEdgePoints_Scatter(subd);
			ccs__CreasedCageVertexPoints_Scatter(subd);

			for (int depth = 1; depth < ccs_MaxDepth(subd); ++depth)
			{
				ccs__FacePoints_Scatter(subd, depth);
				ccs__CreasedEdgePoints_Scatter(subd, depth);
				ccs__CreasedVertexPoints_Scatter(subd, depth);
			}
		}
		else
		{
			ccs__CageFacePoints_Scatter(subd);
			ccs__CageEdgePoints_Scatter(subd);
			ccs__CageVertexPoints_Scatter(subd);

			for (int depth = 1; depth < ccs_MaxDepth(subd); ++depth)
			{
				ccs__FacePoints_Scatter(subd, depth);
				ccs__EdgePoints_Scatter(subd, depth);
				ccs__VertexPoints_Scatter(subd, depth);
			}
		}
	}

	public static void ccs_RefineVertexPoints_Gather(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;

		if (true/*cm_IsCreased(cage)*/)
		{
			ccs__CageFacePoints_Gather(subd);
			ccs__CreasedCageEdgePoints_Gather(subd);
			ccs__CreasedCageVertexPoints_Gather(subd);

			for (int depth = 1; depth < ccs_MaxDepth(subd); ++depth)
			{
				ccs__FacePoints_Gather(subd, depth);
				ccs__CreasedEdgePoints_Gather(subd, depth);
				ccs__CreasedVertexPoints_Gather(subd, depth);
			}
		}
		else
		{
			ccs__CageFacePoints_Gather(subd);
			ccs__CageEdgePoints_Gather(subd);
			ccs__CageVertexPoints_Gather(subd);

			for (int depth = 1; depth < ccs_MaxDepth(subd); ++depth)
			{
				ccs__FacePoints_Gather(subd, depth);
				ccs__EdgePoints_Gather(subd, depth);
				ccs__VertexPoints_Gather(subd, depth);
			}
		}
	}


	/*******************************************************************************
	 * RefineCageHalfedges -- Applies halfedge refinement rules on the cage mesh
	 *
	 * This routine computes the halfedges of the control cage after one subdivision
	 * step and stores them in the subd.
	 *
	 */
	public static void ccs__RefineCageHalfedges(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int vertexCount = ccm_VertexCount(cage);
		int edgeCount = ccm_EdgeCount(cage);
		int faceCount = ccm_FaceCount(cage);
		int halfedgeCount = ccm_HalfedgeCount(cage);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int twinID = ccm_HalfedgeTwinID(cage, halfedgeID);
			int prevID = ccm_HalfedgePrevID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			int faceID = ccm_HalfedgeFaceID(cage, halfedgeID);
			int edgeID = ccm_HalfedgeEdgeID(cage, halfedgeID);
			int prevEdgeID = ccm_HalfedgeEdgeID(cage, prevID);
			int prevTwinID = ccm_HalfedgeTwinID(cage, prevID);
			int vertexID = ccm_HalfedgeVertexID(cage, halfedgeID);
			int twinNextID =
				twinID >= 0 ? ccm_HalfedgeNextID(cage, twinID) : -1;

			// twinIDs
			subd.halfedges[(4 * halfedgeID + 0)].twinID = 4 * twinNextID + 3;
			subd.halfedges[(4 * halfedgeID + 1)].twinID = 4 * nextID + 2;
			subd.halfedges[(4 * halfedgeID + 2)].twinID = 4 * prevID + 1;
			subd.halfedges[(4 * halfedgeID + 3)].twinID = 4 * prevTwinID + 0;

			// edgeIDs
			subd.halfedges[(4 * halfedgeID + 0)].edgeID = 2 * edgeID + (halfedgeID > twinID ? 0 : 1);
			subd.halfedges[(4 * halfedgeID + 1)].edgeID = 2 * edgeCount + halfedgeID;
			subd.halfedges[(4 * halfedgeID + 2)].edgeID = 2 * edgeCount + prevID;
			subd.halfedges[(4 * halfedgeID + 3)].edgeID = 2 * prevEdgeID + (prevID > prevTwinID ? 1 : 0);

			// vertexIDs
			subd.halfedges[(4 * halfedgeID + 0)].vertexID = vertexID;
			subd.halfedges[(4 * halfedgeID + 1)].vertexID = vertexCount + faceCount + edgeID;
			subd.halfedges[(4 * halfedgeID + 2)].vertexID = vertexCount + faceID;
			subd.halfedges[(4 * halfedgeID + 3)].vertexID = vertexCount + faceCount + prevEdgeID;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineHalfedges -- Applies halfedge refinement on the subd
	 *
	 * This routine computes the halfedges of the next subd level.
	 *
	 */
	public static void ccs__RefineHalfedges(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int vertexCount = ccm_VertexCountAtDepth_Fast(cage, depth);
		int edgeCount = ccm_EdgeCountAtDepth_Fast(cage, depth);
		int faceCount = ccm_FaceCountAtDepth_Fast(cage, depth);
		int stride = ccs_CumulativeHalfedgeCountAtDepth(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int twinID = ccs_HalfedgeTwinID(subd, halfedgeID, depth);
			int prevID = ccm_HalfedgePrevID_Quad(halfedgeID);
			int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);
			int faceID = ccm_HalfedgeFaceID_Quad(halfedgeID);
			int edgeID = ccs_HalfedgeEdgeID(subd, halfedgeID, depth);
			int vertexID = ccs_HalfedgeVertexPointID(subd, halfedgeID, depth);
			int prevEdgeID = ccs_HalfedgeEdgeID(subd, prevID, depth);
			int prevTwinID = ccs_HalfedgeTwinID(subd, prevID, depth);
			int twinNextID = ccm_HalfedgeNextID_Quad(twinID);

			// twinIDs
			subd.halfedges[stride + (4 * halfedgeID + 0)].twinID = 4 * twinNextID + 3;
			subd.halfedges[stride + (4 * halfedgeID + 1)].twinID = 4 * nextID + 2;
			subd.halfedges[stride + (4 * halfedgeID + 2)].twinID = 4 * prevID + 1;
			subd.halfedges[stride + (4 * halfedgeID + 3)].twinID = 4 * prevTwinID + 0;

			// edgeIDs
			subd.halfedges[stride + (4 * halfedgeID + 0)].edgeID = 2 * edgeID + (halfedgeID > twinID ? 0 : 1);
			subd.halfedges[stride + (4 * halfedgeID + 1)].edgeID = 2 * edgeCount + halfedgeID;
			subd.halfedges[stride + (4 * halfedgeID + 2)].edgeID = 2 * edgeCount + prevID;
			subd.halfedges[stride + (4 * halfedgeID + 3)].edgeID = 2 * prevEdgeID + (prevID > prevTwinID ? 1 : 0);

			// vertexIDs
			subd.halfedges[stride + (4 * halfedgeID + 0)].vertexID = vertexID;
			subd.halfedges[stride + (4 * halfedgeID + 1)].vertexID = vertexCount + faceCount + edgeID;
			subd.halfedges[stride + (4 * halfedgeID + 2)].vertexID = vertexCount + faceID;
			subd.halfedges[stride + (4 * halfedgeID + 3)].vertexID = vertexCount + faceCount + prevEdgeID;
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineHalfedges
	 *
	 */
	public static void ccs_RefineHalfedges(cc_Subd subd)
	{
		int maxDepth = ccs_MaxDepth(subd);

		ccs__RefineCageHalfedges(subd);

		for (int depth = 1; depth < maxDepth; ++depth)
		{
			ccs__RefineHalfedges(subd, depth);
		}
	}


	/*******************************************************************************
	 * RefineCageVertexUvs -- Refines UVs of the cage mesh
	 *
	 * This routine computes the UVs of the control cage after one subdivision
	 * step and stores them in the subd. Note that since UVs are not linked to
	 * the topology of the mesh, we store the results of the UV computation
	 * within the halfedge buffer.
	 *
	 */
	public static void ccs__RefineCageVertexUvs(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCount(cage);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int prevID = ccm_HalfedgePrevID(cage, halfedgeID);
			int nextID = ccm_HalfedgeNextID(cage, halfedgeID);
			float2 uv = ccm_HalfedgeVertexUv(cage, halfedgeID);
			float2 nextUv = ccm_HalfedgeVertexUv(cage, nextID);
			float2 prevUv = ccm_HalfedgeVertexUv(cage, prevID);
			float2 edgeUv, prevEdgeUv;
			float2 faceUv = uv;
			int m = 1;

			edgeUv = math.lerp(uv, nextUv, 0.5f);
			prevEdgeUv = math.lerp(uv, prevUv, 0.5f);

			for (int halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeID);
						 halfedgeIt != halfedgeID;
						 halfedgeIt = ccm_HalfedgeNextID(cage, halfedgeIt))
			{
				float2 uv2 = ccm_HalfedgeVertexUv(cage, halfedgeIt);
				faceUv += uv2;
				++m;
			}
			faceUv /= (float)m;

			subd.halfedges[(4 * halfedgeID + 0)].uvID = cc__EncodeUv(uv);
			subd.halfedges[(4 * halfedgeID + 1)].uvID = cc__EncodeUv(edgeUv);
			subd.halfedges[(4 * halfedgeID + 2)].uvID = cc__EncodeUv(faceUv);
			subd.halfedges[(4 * halfedgeID + 3)].uvID = cc__EncodeUv(prevEdgeUv);
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineVertexUvs -- Applies UV refinement on the subd
	 *
	 * This routine computes the UVs of the next subd level.
	 *
	 */
	public static void ccs__RefineVertexUvs(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
		int stride = ccs_CumulativeHalfedgeCountAtDepth(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int prevID = ccm_HalfedgePrevID_Quad(halfedgeID);
			int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);
			float2 uv = ccs_HalfedgeVertexUv(subd, halfedgeID, depth);
			float2 nextUv = ccs_HalfedgeVertexUv(subd, nextID, depth);
			float2 prevUv = ccs_HalfedgeVertexUv(subd, prevID, depth);
			float2 edgeUv, prevEdgeUv;
			float2 faceUv = uv;

			edgeUv = math.lerp(uv, nextUv, 0.5f);
			prevEdgeUv = math.lerp(uv, prevUv, 0.5f);

			for (int halfedgeIt = ccs_HalfedgeNextID(subd, halfedgeID, depth);
						 halfedgeIt != halfedgeID;
						 halfedgeIt = ccs_HalfedgeNextID(subd, halfedgeIt, depth))
			{
				float2 uv2 = ccs_HalfedgeVertexUv(subd, halfedgeIt, depth);
				faceUv += uv2;
			}
			faceUv /= 4.0f;

			subd.halfedges[stride + (4 * halfedgeID + 0)].uvID = ccs__HalfedgeVertexUvID(subd, halfedgeID, depth);
			subd.halfedges[stride + (4 * halfedgeID + 1)].uvID = cc__EncodeUv(edgeUv);
			subd.halfedges[stride + (4 * halfedgeID + 2)].uvID = cc__EncodeUv(faceUv);
			subd.halfedges[stride + (4 * halfedgeID + 3)].uvID = cc__EncodeUv(prevEdgeUv);
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineUvs
	 *
	 */
	public static void ccs_RefineVertexUvs(cc_Subd subd)
	{
		if (ccm_UvCount(subd.cage) > 0)
		{
			int maxDepth = ccs_MaxDepth(subd);

			ccs__RefineCageVertexUvs(subd);

			for (int depth = 1; depth < maxDepth; ++depth)
			{
				ccs__RefineVertexUvs(subd, depth);
			}
		}
	}


	/*******************************************************************************
	 * RefineCageCreases -- Applies crease subdivision on the cage mesh
	 *
	 * This routine computes the creases of the control cage after one subdivision
	 * step and stores them in the subd.
	 *
	 */
	public static void ccs__RefineCageCreases(cc_Subd subd)
	{
		cc_Mesh cage = subd.cage;
		int edgeCount = ccm_EdgeCount(cage);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int nextID = ccm_CreaseNextID(cage, edgeID);
			int prevID = ccm_CreasePrevID(cage, edgeID);
			bool t1 = ccm_CreasePrevID(cage, nextID) == edgeID;
			bool t2 = ccm_CreaseNextID(cage, prevID) == edgeID;
			float thisS = 3.0f * ccm_CreaseSharpness(cage, edgeID);
			float nextS = ccm_CreaseSharpness(cage, nextID);
			float prevS = ccm_CreaseSharpness(cage, prevID);

			// next rule
			subd.creases[(2 * edgeID + 0)].nextID = 2 * edgeID + 1;
			subd.creases[(2 * edgeID + 1)].nextID = 2 * nextID + (t1 ? 0 : 1);

			// prev rule
			subd.creases[(2 * edgeID + 0)].prevID = 2 * prevID + (t2 ? 1 : 0);
			subd.creases[(2 * edgeID + 1)].prevID = 2 * edgeID + 0;

			// sharpness rule
			subd.creases[(2 * edgeID + 0)].sharpness = math.max(0.0f, (prevS + thisS) / 4.0f - 1.0f);
			subd.creases[(2 * edgeID + 1)].sharpness = math.max(0.0f, (thisS + nextS) / 4.0f - 1.0f);
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineCreases -- Applies crease subdivision on the subd
	 *
	 * This routine computes the topology of the next subd level.
	 *
	 */
	public static void ccs__RefineCreases(cc_Subd subd, int depth)
	{
		cc_Mesh cage = subd.cage;
		int creaseCount = ccm_CreaseCountAtDepth(cage, depth);
		int stride = ccs_CumulativeCreaseCountAtDepth(cage, depth);

		// CC_PARALLEL_FOR // TODO: TRY TO USE JOB SYSTEM

		for (int edgeID = 0; edgeID < creaseCount; ++edgeID)
		{
			int nextID = ccs_CreaseNextID_Fast(subd, edgeID, depth);
			int prevID = ccs_CreasePrevID_Fast(subd, edgeID, depth);
			bool t1 = ccs_CreasePrevID_Fast(subd, nextID, depth) == edgeID;
			bool t2 = ccs_CreaseNextID_Fast(subd, prevID, depth) == edgeID;
			float thisS = 3.0f * ccs_CreaseSharpness_Fast(subd, edgeID, depth);
			float nextS = ccs_CreaseSharpness_Fast(subd, nextID, depth);
			float prevS = ccs_CreaseSharpness_Fast(subd, prevID, depth);

			// next rule
			subd.creases[stride + (2 * edgeID + 0)].nextID = 2 * edgeID + 1;
			subd.creases[stride + (2 * edgeID + 1)].nextID = 2 * nextID + (t1 ? 0 : 1);

			// prev rule
			subd.creases[stride + (2 * edgeID + 0)].prevID = 2 * prevID + (t2 ? 1 : 0);
			subd.creases[stride + (2 * edgeID + 1)].prevID = 2 * edgeID + 0;

			// sharpness rule
			subd.creases[stride + (2 * edgeID + 0)].sharpness = math.max(0.0f, (prevS + thisS) / 4.0f - 1.0f);
			subd.creases[stride + (2 * edgeID + 1)].sharpness = math.max(0.0f, (thisS + nextS) / 4.0f - 1.0f);
		}
		//CC_BARRIER
	}


	/*******************************************************************************
	 * RefineCreases
	 *
	 */
	public static void ccs_RefineCreases(cc_Subd subd)
	{
		int maxDepth = ccs_MaxDepth(subd);

		ccs__RefineCageCreases(subd);

		for (int depth = 1; depth < maxDepth; ++depth)
		{
			ccs__RefineCreases(subd, depth);
		}
	}


	/*******************************************************************************
	 * Refine -- Computes and stores the result of Catmull Clark subdivision.
	 *
	 * The subdivision is computed down to the maxDepth parameter.
	 *
	 */
	public static void ccs__RefineTopology(cc_Subd subd)
	{
		ccs_RefineHalfedges(subd);
		ccs_RefineVertexUvs(subd);
		ccs_RefineCreases(subd);
	}

	public static void ccs_Refine_Scatter(cc_Subd subd)
	{
		ccs__RefineTopology(subd);
		ccs_RefineVertexPoints_Scatter(subd);
	}

	public static void ccs_Refine_Gather(cc_Subd subd)
	{
		ccs__RefineTopology(subd);
		ccs_RefineVertexPoints_Gather(subd);
	}


	/*******************************************************************************
	 * Magic -- Generates the magic identifier
	 *
	 * Each cc_Mesh file starts with 8 Bytes that allow us to check if the file
	 * under reading is actually a cc_Mesh file.
	 *
	 */
	public static long ccm__Magic()
	{
		char[] chars = new char[] { 'c', 'c', '_', 'M', 'e', 's', 'h', '1' };
		byte[] bytes = new byte[8];
		for (int i = 0; i < bytes.Length; i++)
			bytes[i] = System.Convert.ToByte(chars[i]);

		return System.BitConverter.ToInt64(bytes, 0);
	}


	/*******************************************************************************
	 * Header File Data Structure
	 *
	 * This represents the header we use to uniquely identify the cc_Mesh files
	 * and provide the fundamental information to properly decode the rest of the
	 * file.
	 *
	 */
	public struct ccm__Header
	{
		public long magic;
		public int vertexCount;
		public int uvCount;
		public int halfedgeCount;
		public int edgeCount;
		public int faceCount;
	}


	/*******************************************************************************
	 * ReadHeader -- Reads a tt_Texture file header from an input stream
	 *
	 */
	public static bool ccm__ReadHeader(BinaryReader stream, out ccm__Header header)
	{
		header = new ccm__Header();
		try
		{
			header.magic = stream.ReadInt64();
			header.vertexCount = stream.ReadInt32();
			header.uvCount = stream.ReadInt32();
			header.halfedgeCount = stream.ReadInt32();
			header.edgeCount = stream.ReadInt32();
			header.faceCount = stream.ReadInt32();
		}
		catch
		{
			Debug.Log("cc: ccm__ReadHeader failed");
			return false;
		}
		return header.magic == ccm__Magic();
	}


	/*******************************************************************************
	 * ReadData -- Loads mesh data
	 *
	 */
	public static bool ccm__ReadData(cc_Mesh mesh, BinaryReader stream)
	{
		int vertexCount = ccm_VertexCount(mesh);
		int uvCount = ccm_UvCount(mesh);
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		int creaseCount = ccm_CreaseCount(mesh);
		int edgeCount = ccm_EdgeCount(mesh);
		int faceCount = ccm_FaceCount(mesh);

		try
		{
			stream.ReadInt32(); // TODO: FIND OUT WHY NECESSARY
			for (int i = 0; i < vertexCount; i++)
				mesh.vertexToHalfedgeIDs[i] = stream.ReadInt32();
			for (int i = 0; i < edgeCount; i++)
				mesh.edgeToHalfedgeIDs[i] = stream.ReadInt32();
			for (int i = 0; i < faceCount; i++)
				mesh.faceToHalfedgeIDs[i] = stream.ReadInt32();
			for (int i = 0; i < vertexCount; i++)
				mesh.vertexPoints[i] = new float3(stream.ReadSingle(), stream.ReadSingle(), stream.ReadSingle());
			for (int i = 0; i < uvCount; i++)
				mesh.uvs[i] = new float2(stream.ReadSingle(), stream.ReadSingle());
			for (int i = 0; i < creaseCount; i++)
			{
				mesh.creases[i].nextID = stream.ReadInt32();
				mesh.creases[i].prevID = stream.ReadInt32();
				mesh.creases[i].sharpness = stream.ReadSingle();
			}
			for (int i = 0; i < halfedgeCount; i++)
			{
				mesh.halfedges[i].twinID = stream.ReadInt32();
				mesh.halfedges[i].nextID = stream.ReadInt32();
				mesh.halfedges[i].prevID = stream.ReadInt32();
				mesh.halfedges[i].faceID = stream.ReadInt32();
				mesh.halfedges[i].edgeID = stream.ReadInt32();
				mesh.halfedges[i].vertexID = stream.ReadInt32();
				mesh.halfedges[i].uvID = stream.ReadInt32();
			}
		}
		catch
		{
			return false;
		}
		return true;
	}


	/*******************************************************************************
	 * Load -- Loads a mesh from a CCM file
	 *
	 */
	public static cc_Mesh ccm_Load(string filename)
	{
		cc_Mesh mesh = new cc_Mesh();

		if (File.Exists(filename) == false)
		{
			Debug.Log("cc: file open failed");
			return null;
		}

		BinaryReader stream = new BinaryReader(File.Open(filename, FileMode.Open));

		ccm__Header header;
		if (!ccm__ReadHeader(stream, out header))
		{
			Debug.Log("cc: unsupported file");
			stream.Close();
			return null;
		}

		mesh = ccm_Create(header.vertexCount,
						  header.uvCount,
						  header.halfedgeCount,
						  header.edgeCount,
						  header.faceCount);
		if (!ccm__ReadData(mesh, stream))
		{
			Debug.Log("cc: data reading failed");
			stream.Close();
			return null;
		}
		stream.Close();

		if (header.uvCount == 0)
		{
			header.uvCount = header.vertexCount;
			mesh.uvCount = mesh.vertexCount;
			mesh.uvs = new float2[mesh.uvCount];
		}

		return mesh;
	}

	/*******************************************************************************
	 * Save -- Save a mesh to a file
	 *
	 */
	public static bool ccm_Save(cc_Mesh mesh, string filename)
	{
		int vertexCount = ccm_VertexCount(mesh);
		int uvCount = ccm_UvCount(mesh);
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		int creaseCount = ccm_CreaseCount(mesh);
		int edgeCount = ccm_EdgeCount(mesh);
		int faceCount = ccm_FaceCount(mesh);

		ccm__Header header;
		header.magic = ccm__Magic();
		header.vertexCount = ccm_VertexCount(mesh);
		header.uvCount = ccm_UvCount(mesh);
		header.halfedgeCount = ccm_HalfedgeCount(mesh);
		header.edgeCount = ccm_EdgeCount(mesh);
		header.faceCount = ccm_FaceCount(mesh);

		BinaryWriter stream = new BinaryWriter(File.Open(filename, FileMode.Create));

		// Write header
		try
		{
			stream.Write(header.magic);
			stream.Write(header.vertexCount);
			stream.Write(header.uvCount);
			stream.Write(header.halfedgeCount);
			stream.Write(header.edgeCount);
			stream.Write(header.faceCount);
		}
		catch
		{
			Debug.Log("cc: ccm__WriteHeader failed");
			return false;
		}

		// Write data
		try
		{
			stream.Write(0); // TODO: FIND OUT WHY NECESSARY
			for (int i = 0; i < vertexCount; i++)
				stream.Write(mesh.vertexToHalfedgeIDs[i]);
			for (int i = 0; i < edgeCount; i++)
				stream.Write(mesh.edgeToHalfedgeIDs[i]);
			for (int i = 0; i < faceCount; i++)
				stream.Write(mesh.faceToHalfedgeIDs[i]);
			for (int i = 0; i < vertexCount; i++)
			{
				stream.Write(mesh.vertexPoints[i].x); stream.Write(mesh.vertexPoints[i].y); stream.Write(mesh.vertexPoints[i].z);
			}
			for (int i = 0; i < uvCount; i++)
			{
				stream.Write(mesh.uvs[i].x); stream.Write(mesh.uvs[i].y);
			}
			for (int i = 0; i < creaseCount; i++)
			{
				stream.Write(mesh.creases[i].nextID);
				stream.Write(mesh.creases[i].prevID);
				stream.Write(mesh.creases[i].sharpness);
			}
			for (int i = 0; i < halfedgeCount; i++)
			{
				stream.Write(mesh.halfedges[i].twinID);
				stream.Write(mesh.halfedges[i].nextID);
				stream.Write(mesh.halfedges[i].prevID);
				stream.Write(mesh.halfedges[i].faceID);
				stream.Write(mesh.halfedges[i].edgeID);
				stream.Write(mesh.halfedges[i].vertexID);
				stream.Write(mesh.halfedges[i].uvID);
			}
		}
		catch
		{
			Debug.Log("cc: data writing failed");
			return false;
		}

		stream.Close();
		return true;
	}

	/*******************************************************************************
	 * ExportToObj -- Exports subd to the OBJ file format
	 *
	 */
	public static void ExportToObj(cc_Subd subd, int depth, string filename)
	{
		cc_Mesh cage = subd.cage;
		int vertexPointCount = ccm_VertexCountAtDepth(cage, depth);
		int faceCount = ccm_FaceCountAtDepth(cage, depth);


		StreamWriter stream = new StreamWriter(File.Open(filename, FileMode.Create));

		// write vertices
		stream.Write("# Vertices\n");
		if (depth == 0)
		{
			int vertexUvCount = ccm_UvCount(cage);
			for (int vertexID = 0; vertexID < vertexPointCount; ++vertexID)
			{
				float3 v = ccm_VertexPoint(cage, vertexID);
				stream.Write("v {0} {1} {2}\n", v[0], v[1], v[2]);
			}
			for (int vertexID = 0; vertexID < vertexUvCount; ++vertexID)
			{
				float2 v = ccm_Uv(cage, vertexID);
				stream.Write("vt {0} {1}\n", v[0], v[1]);
			}
		}
		else
		{
			int halfedgeCount = ccm_HalfedgeCountAtDepth(cage, depth);
			for (int vertexID = 0; vertexID < vertexPointCount; ++vertexID)
			{
				float3 v = ccs_VertexPoint(subd, vertexID);
				stream.Write("v {0} {1} {2}\n", v[0], v[1], v[2]);
			}
			for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
			{
				float2 uv = ccs_HalfedgeVertexUv(subd, halfedgeID, depth);
				stream.Write("vt {0} {1}\n", uv[0], uv[1]);
			}
		}
		stream.Write("\n");
		// write topology
		stream.Write("# Topology\n");
		if (depth == 0)
		{
			for (int faceID = 0; faceID < faceCount; ++faceID)
			{
				int halfEdgeID = ccm_FaceToHalfedgeID(cage, faceID);
				stream.Write(
					"f {0}/{1}",
					ccm_HalfedgeVertexID(cage, halfEdgeID) + 1,
					ccm_HalfedgeUvID(cage, halfEdgeID) + 1);
				for (int halfEdgeIt = ccm_HalfedgeNextID(cage, halfEdgeID);
					halfEdgeIt != halfEdgeID;
					halfEdgeIt = ccm_HalfedgeNextID(cage, halfEdgeIt))
				{
					stream.Write(
						" {0}/{1}",
						ccm_HalfedgeVertexID(cage, halfEdgeIt) + 1,
						ccm_HalfedgeUvID(cage, halfEdgeIt) + 1);
				}
				stream.Write("\n");
			}
		}
		else
		{
			for (int faceID = 0; faceID < faceCount; ++faceID)
			{
				stream.Write(
					"f {0}/{1} {2}/{3} {4}/{5} {6}/{7}\n",
					ccs_HalfedgeVertexPointID(subd, 4 * faceID + 0, depth) + 1,
					4 * faceID + 1,
					ccs_HalfedgeVertexPointID(subd, 4 * faceID + 1, depth) + 1,
					4 * faceID + 2,
					ccs_HalfedgeVertexPointID(subd, 4 * faceID + 2, depth) + 1,
					4 * faceID + 3,
					ccs_HalfedgeVertexPointID(subd, 4 * faceID + 3, depth) + 1,
					4 * faceID + 4);
			}
			stream.Write("\n");
		}
		stream.Close();
	}






	/*******************************************************************************
	 * ComputeTwins -- Computes the twin of each half edge
	 *
	 * This routine is what effectively converts a traditional "indexed mesh"
	 * into a halfedge mesh (in the case where all the primitives are the same).
	 *
	 */
	public struct TwinComputationData
	{
		public int halfedgeID;
		public ulong hashID;
	}

	public static int BinarySearch(
		TwinComputationData[] data,
		ulong hashID,
		int beginID,
		int endID)
	{
		int midID;

		if (beginID > endID)
			return -1; // not found

		midID = (beginID + endID) / 2;

		if (data[midID].hashID == hashID)
		{
			return data[midID].halfedgeID;
		}
		else if (hashID > data[midID].hashID)
		{
			return BinarySearch(data, hashID, midID + 1, endID);
		}
		else
		{
			return BinarySearch(data, hashID, beginID, midID - 1);
		}
	}

	public static TwinComputationData[] ComputeTwins(cc_Mesh mesh)
	{
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		ulong vertexCount = (ulong)ccm_VertexCount(mesh);
		Debug.Assert(vertexCount < 0xFFFFFFFFu, "CatmullClarkMeshRenderer: Unsupported Geometry, too many vertices");
		TwinComputationData[] table = new TwinComputationData[halfedgeCount];

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
			ulong v0 = (ulong)ccm_HalfedgeVertexID(mesh, halfedgeID);
			ulong v1 = (ulong)ccm_HalfedgeVertexID(mesh, nextID);

			table[halfedgeID].halfedgeID = halfedgeID;
			table[halfedgeID].hashID = v0 + vertexCount * v1;
		}

		System.Array.Sort(table, (a, b) => a.hashID.CompareTo(b.hashID));

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
			ulong v0 = (ulong)ccm_HalfedgeVertexID(mesh, halfedgeID);
			ulong v1 = (ulong)ccm_HalfedgeVertexID(mesh, nextID);
			ulong twinHashID = v1 + vertexCount * v0;
			int twinID = BinarySearch(table, twinHashID, 0, halfedgeCount - 1);

			mesh.halfedges[halfedgeID].twinID = twinID;
		}

		return table;
	}


	/*******************************************************************************
	 * LoadFaceMappings -- Computes the mappings for the faces of the mesh
	 *
	 */
	public static int FaceScroll(int id, int direction, int maxValue)
	{
		int n = maxValue - 1;
		int d = direction;
		int u = (d + 1) >> 1; // in [0, 1]
		int un = u * n; // precomputation

		return (id == un) ? (n - un) : (id + d);
	}

	public static int ScrollFaceHalfedgeID(
		int halfedgeID,
		int halfedgeFaceBeginID,
		int halfedgeFaceEndID,
		int direction)
	{
		int faceHalfedgeCount = halfedgeFaceEndID - halfedgeFaceBeginID;
		int localHalfedgeID = halfedgeID - halfedgeFaceBeginID;
		int nextHalfedgeID = FaceScroll(localHalfedgeID,
										direction,
										faceHalfedgeCount);

		return halfedgeFaceBeginID + nextHalfedgeID;
	}

	public static void LoadFaceMappings(cc_Mesh mesh, List<MeshTopology> submeshTopologies, List<int[]> submeshIndices)
	{
		// Initialize bitfield for face mappings
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		BitField faceIterator = new BitField(halfedgeCount + 1);
		int halfedgeOffset = 0;
		for (int submeshID = 0; submeshID < submeshIndices.Count; submeshID++)
		{
			if (submeshTopologies[submeshID] == MeshTopology.Triangles)
				for (int triangleID = 0; triangleID < submeshIndices[submeshID].Length / 3; triangleID++)
					faceIterator.SetBit(triangleID * 3 + halfedgeOffset, 1);
			else
				for (int quadID = 0; quadID < submeshIndices[submeshID].Length / 4; quadID++)
					faceIterator.SetBit(quadID * 4 + halfedgeOffset, 1);
			halfedgeOffset += submeshIndices[submeshID].Length;
		}
		faceIterator.SetBit(halfedgeCount, 1);
		faceIterator.Reduce();

		int faceCount = faceIterator.BitCount() - 1;

		mesh.faceToHalfedgeIDs = new int[faceCount];
		mesh.faceCount = faceCount;

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int tmp = faceIterator.EncodeBit(halfedgeID);
			int faceID = tmp - (faceIterator.GetBit(halfedgeID) ^ 1);

			mesh.halfedges[halfedgeID].faceID = faceID;
		}

		for (int faceID = 0; faceID < faceCount; ++faceID)
		{
			mesh.faceToHalfedgeIDs[faceID] = faceIterator.DecodeBit(faceID);
		}

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int faceID = mesh.halfedges[halfedgeID].faceID;
			int beginID = faceIterator.DecodeBit(faceID);
			int endID = faceIterator.DecodeBit(faceID + 1);
			int nextID = ScrollFaceHalfedgeID(halfedgeID, beginID, endID, +1);
			int prevID = ScrollFaceHalfedgeID(halfedgeID, beginID, endID, -1);

			mesh.halfedges[halfedgeID].nextID = nextID;
			mesh.halfedges[halfedgeID].prevID = prevID;
		}
	}

	/*******************************************************************************
	 * LoadEdgeMappings -- Computes the mappings for the edges of the mesh
	 *
	 * Catmull-Clark subdivision requires access to the edges of an input mesh.
	 * Since we are dealing with a half-edge representation, we virtually
	 * have to iterate the half-edges in a sparse way (an edge is a pair of
	 * neighboring half-edges in the general case, except for boundary edges
	 * where it only consists of a single half-edge).
	 * This function builds a data-structure that allows to do just that:
	 * for each halfedge pair, we only consider the one that has the largest
	 * halfedgeID. This allows to treat boundary and regular edges seamlessly.
	 *
	 */
	public static void LoadEdgeMappings(cc_Mesh mesh)
	{
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		BitField edgeIterator = new BitField(halfedgeCount);
		int edgeCount;

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);
			int bitValue = halfedgeID > twinID ? 1 : 0;

			edgeIterator.SetBit(halfedgeID, bitValue);
		}

		edgeIterator.Reduce();
		edgeCount = edgeIterator.BitCount();

		mesh.edgeToHalfedgeIDs = new int[edgeCount];
		mesh.edgeCount = edgeCount;

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);
			int bitID = halfedgeID > twinID ? halfedgeID : twinID;

			mesh.halfedges[halfedgeID].edgeID = edgeIterator.EncodeBit(bitID);
		}

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			mesh.edgeToHalfedgeIDs[edgeID] = edgeIterator.DecodeBit(edgeID);
		}
	}

	/*******************************************************************************
	 * LoadVertexHalfedges -- Computes an iterator over one half-edge per vertex
	 *
	 * Catmull-Clark subdivision requires access to the half-edges that surround
	 * the vertices of an input mesh.
	 * This function determines a half-edge ID that starts from a
	 * given vertex within that vertex. We distinguish two cases:
	 * 1- If the vertex is a lying on a boundary, we stored the halfedge that
	 * allows for iteration in the forward sense.
	 * 2- Otherwise we store the largest half-edge ID.
	 *
	 */
	public static void LoadVertexHalfedges(cc_Mesh mesh, out bool crashed)
	{
		crashed = false;
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		int vertexCount = ccm_VertexCount(mesh);

		mesh.vertexToHalfedgeIDs = new int[vertexCount];

		for (int halfedgeID = 0; halfedgeID < halfedgeCount; ++halfedgeID)
		{
			int vertexID = ccm_HalfedgeVertexID(mesh, halfedgeID);
			int maxHalfedgeID = halfedgeID;
			int boundaryHalfedgeID = halfedgeID;
			int iterator;

			int crashCounter = 0;
			for (iterator = ccm_NextVertexHalfedgeID(mesh, halfedgeID);
				 iterator >= 0 && iterator != halfedgeID && crashCounter < 1000;
				 iterator = ccm_NextVertexHalfedgeID(mesh, iterator))
			{
				maxHalfedgeID = maxHalfedgeID > iterator ? maxHalfedgeID : iterator;
				boundaryHalfedgeID = iterator;
				crashCounter++;
			}
			if (crashCounter >= 1000)
			{
				crashed = true;
				Debug.LogError("CatmullClarkMeshRenderer: Unsupported Geometry");
				return;
			}

			// affect max half-edge ID to vertex
			if /*boundary involved*/ (iterator < 0)
			{
				if (halfedgeID == boundaryHalfedgeID)
				{
					mesh.vertexToHalfedgeIDs[vertexID] = boundaryHalfedgeID;
				}
			}
			else
			{
				if (halfedgeID == maxHalfedgeID)
				{
					mesh.vertexToHalfedgeIDs[vertexID] = maxHalfedgeID;
				}
			}
		}
	}

	/*******************************************************************************
	 * MakeBoundariesSharp -- Tags boundary edges as sharp
	 *
	 * Following the Pixar standard, we tag boundary halfedges as sharp.
	 * See "Subdivision Surfaces in Character Animation" by DeRose et al.
	 * Note that we tag the sharpness value to 16 as subdivision can't go deeper
	 * without overflowing 32-bit integers.
	 *
	 */
	public static void InitCreasesAndMakeBoundariesSharp(cc_Mesh mesh)
	{
		int edgeCount = ccm_EdgeCount(mesh);

		mesh.creases = new cc_Crease[edgeCount];
		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			mesh.creases[edgeID].nextID = edgeID;
			mesh.creases[edgeID].prevID = edgeID;
		}

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			int halfedgeID = ccm_EdgeToHalfedgeID(mesh, edgeID);
			int twinID = ccm_HalfedgeTwinID(mesh, halfedgeID);

			if (twinID < 0)
			{
				mesh.creases[edgeID].sharpness = 16.0f;
			}
		}
	}

	/*******************************************************************************
	 * ComputeCreaseNeighbors -- Computes the neighbors of each crease
	 *
	 */
	public static void ComputeCreaseNeighbors(cc_Mesh mesh)
	{
		int edgeCount = ccm_EdgeCount(mesh);

		for (int edgeID = 0; edgeID < edgeCount; ++edgeID)
		{
			float sharpness = ccm_CreaseSharpness(mesh, edgeID);

			if (sharpness > 0.0f)
			{
				int halfedgeID = ccm_EdgeToHalfedgeID(mesh, edgeID);
				int nextID = ccm_HalfedgeNextID(mesh, halfedgeID);
				int prevCreaseCount = 0;
				int prevCreaseID = -1;
				int nextCreaseCount = 0;
				int nextCreaseID = -1;
				int halfedgeIt;

				for (halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeID);
					 halfedgeIt != halfedgeID && halfedgeIt >= 0;
					 halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeIt))
				{
					float s = ccm_HalfedgeSharpness(mesh, halfedgeIt);

					if (s > 0.0f)
					{
						prevCreaseID = ccm_HalfedgeEdgeID(mesh, halfedgeIt);
						++prevCreaseCount;
					}
				}

				if (prevCreaseCount == 1 && halfedgeIt == halfedgeID)
				{
					mesh.creases[edgeID].prevID = prevCreaseID;
				}

				if (ccm_HalfedgeSharpness(mesh, nextID) > 0.0f)
				{
					nextCreaseID = ccm_HalfedgeEdgeID(mesh, nextID);
					++nextCreaseCount;
				}

				for (halfedgeIt = ccm_NextVertexHalfedgeID(mesh, nextID);
					 halfedgeIt != nextID && halfedgeIt >= 0;
					 halfedgeIt = ccm_NextVertexHalfedgeID(mesh, halfedgeIt))
				{
					float s = ccm_HalfedgeSharpness(mesh, halfedgeIt);
					int twinID = ccm_HalfedgeTwinID(mesh, halfedgeIt);

					// twin check is to avoid counting for halfedgeID
					if (s > 0.0f && twinID != halfedgeID)
					{
						nextCreaseID = ccm_HalfedgeEdgeID(mesh, halfedgeIt);
						++nextCreaseCount;
					}
				}

				if (nextCreaseCount == 1 && halfedgeIt == nextID)
				{
					mesh.creases[edgeID].nextID = nextCreaseID;
				}
			}
		}
	}

	/*******************************************************************************
	 * LoadFromUnity -- Loads a CCM mesh from a Unity mesh
	 *
	 * Supports both triangle and quad submeshes (see Keep Quads option). Because
	 * polygons with more than 4 vertices are triangulated by the Unity mesh importer,
	 * the subdivision surface won't produce exactly the same result for these polygons.
	 * Currently all submeshes within the mesh are merged into a single halfedge mesh for
	 * simplicity, which means having a different material per submesh is not supported.
	 */
	public static cc_Mesh ccm_LoadFromUnity(Mesh unityMesh, ref int[] unityVertexBufferToCCMWeldedBuffer, bool weldVertices = true)
	{
		// Load Unity Mesh Data
		Vector3[] vertices = vertices = unityMesh.vertices;
		int vertexCount = vertices.Length;
		Vector2[] texcoords0 = unityMesh.uv;
		int uvCount = texcoords0.Length;

		// Load all triangle and quad submeshes' indices
		int halfedgeCount = 0;
		List<int[]> submeshIndices = new List<int[]>();
		List<int[]> submeshIndicesUnwelded = new List<int[]>();
		List<MeshTopology> submeshTopologies = new List<MeshTopology>();
		for (int submeshID = 0; submeshID < unityMesh.subMeshCount; submeshID++)
		{
			MeshTopology topology = unityMesh.GetTopology(submeshID);
			if (topology != MeshTopology.Triangles && topology != MeshTopology.Quads)
				continue;

			// TEMP TEST
			submeshTopologies.Add(topology);
			submeshIndices.Add(unityMesh.GetIndices(submeshID));
			submeshIndicesUnwelded.Add(unityMesh.GetIndices(submeshID));
			halfedgeCount += submeshIndices[submeshIndices.Count - 1].Length;
		}

		// Remove duplicate vertices and create remapping array to keep separate UVs
		// CAUTION : To use only for vertices that were duplicated by Unity for having separate UVs at uv map discontinuities !
		int[] weldMapping = new int[vertexCount];
		if (weldVertices == true)
		{
			BitField vertexIterator = new BitField(vertexCount);
			List<Vector3> weldedVertices = new List<Vector3>();
			for (int i = 0; i < vertexCount; i++)
			{
				bool vertexExists = false;
				for (int j = 0; j < weldedVertices.Count; j++)
				{
					if (weldVertices == true && weldedVertices[j].x == vertices[i].x && weldedVertices[j].y == vertices[i].y && weldedVertices[j].z == vertices[i].z)
					{
						weldMapping[i] = j;
						vertexExists = true;
						break;
					}
				}
				if (vertexExists == false)
				{
					weldMapping[i] = weldedVertices.Count;
					weldedVertices.Add(vertices[i]);
					vertexIterator.SetBit(i, 1);
				}
			}
			vertices = weldedVertices.ToArray();
			vertexCount = vertices.Length;

			// Bake reverse mapping for later usage
			vertexIterator.Reduce();
			unityVertexBufferToCCMWeldedBuffer = new int[vertexCount];
			for (int bitID = 0; bitID < vertexIterator.BitCount(); bitID++)
			{
				int vertexID = vertexIterator.DecodeBit(bitID);
				unityVertexBufferToCCMWeldedBuffer[bitID] = vertexID;
			}

			// Remap indices to account for vertex welding
			for (int submeshID = 0; submeshID < submeshIndices.Count; submeshID++)
				for (int halfedgeID = 0; halfedgeID < submeshIndices[submeshID].Length; halfedgeID++)
					submeshIndices[submeshID][halfedgeID] = weldMapping[submeshIndices[submeshID][halfedgeID]];
		}
		else
		{
			unityVertexBufferToCCMWeldedBuffer = new int[vertexCount];
			for (int i = 0; i < vertexCount; i++)
				unityVertexBufferToCCMWeldedBuffer[i] = i;
		}

		// Initialize Catmull Clark Mesh structure
		cc_Mesh mesh = new cc_Mesh();
		mesh.vertexCount = vertexCount;
		mesh.uvCount = uvCount;
		mesh.halfedgeCount = halfedgeCount;
		mesh.vertexPoints = new float3[vertexCount];
		mesh.uvs = new float2[uvCount];
		mesh.halfedges = new cc_Halfedge[halfedgeCount];

		// Load vertices
		for (int vertexID = 0; vertexID < vertexCount; vertexID++)
			mesh.vertexPoints[vertexID] = vertices[vertexID];

		// Load uvs
		for (int uvID = 0; uvID < texcoords0.Length; uvID++)
			mesh.uvs[uvID] = texcoords0[uvID];

		// Load halfedge indices (merging all submeshes, triangles and quads, into one halfedge mesh) + handle welded vertices with different UVs
		int halfedgeOffset = 0;
		for (int submeshID = 0; submeshID < submeshIndices.Count; submeshID++)
		{
			for (int halfedgeID = 0; halfedgeID < submeshIndices[submeshID].Length; halfedgeID++)
			{
				mesh.halfedges[halfedgeID + halfedgeOffset].vertexID = submeshIndices[submeshID][halfedgeID];
				mesh.halfedges[halfedgeID + halfedgeOffset].uvID = submeshIndicesUnwelded[submeshID][halfedgeID];
			}
			halfedgeOffset += submeshIndices[submeshID].Length;
		}

		LoadFaceMappings(mesh, submeshTopologies, submeshIndices);
		ComputeTwins(mesh);
		LoadEdgeMappings(mesh);

		bool crashed;
		LoadVertexHalfedges(mesh, out crashed);
		if (crashed == true)
			return null;

		InitCreasesAndMakeBoundariesSharp(mesh);
		// Load creases here
		ComputeCreaseNeighbors(mesh);

		return mesh;
	}

	/*******************************************************************************
	 * BitField class, helper for ccm_LoadFromUnity
	 * 
	 */
	public class BitField
	{
		public int[] m_buffer;


		/*******************************************************************************
		 * Bitfield Constructor
		 *
		 */
		public BitField(int size)
		{
			int sizePowerOfTwo = NextPowerOfTwo(size);
			m_buffer = new int[2 * sizePowerOfTwo];
			m_buffer[0] = sizePowerOfTwo;
		}

		/*******************************************************************************
		 * Size -- Returns the size of the bitfield
		 *
		 */
		public int Size()
		{
			return m_buffer[0];
		}

		/*******************************************************************************
		 * NextPowerOfTwo -- Returns the upper power of two value
		 *
		 * if the input is already a power of two, its value is returned.
		 *
		 */
		private static int NextPowerOfTwo(int x)
		{
			x--;
			x |= x >> 1;
			x |= x >> 2;
			x |= x >> 4;
			x |= x >> 8;
			x |= x >> 16;
			x++;

			return x;
		}

		/*******************************************************************************
		 * FindMSB -- Returns the position of the most significant bit
		 *
		 */
		static private int FindMSB(int x)
		{
			int msb = 0;

			while (x > 1u)
			{
				++msb;
				x = x >> 1;
			}

			return msb;
		}

		/*******************************************************************************
		 * SetBit -- Set a specific bit to either 0 or 1 in the bitfield
		 *
		 */
		public void SetBit(int bitID, int bitValue)
		{
			int offset = Size();
			m_buffer[offset + bitID] = bitValue;
		}

		/*******************************************************************************
		 * GetBit -- Get a specific bit in the bitfield
		 *
		 */
		public int GetBit(int bitID)
		{
			int offset = Size();
			return m_buffer[offset + bitID];
		}

		/*******************************************************************************
		 * BitCount -- Returns the number of bits set to one in the bit field
		 *
		 */
		public int BitCount()
		{
			return m_buffer[1];
		}

		/*******************************************************************************
		 * DecodeNode -- Returns the leaf node associated to index nodeID
		 *
		 * This is procedure is for iterating over the one-valued bits.
		 *
		 */
		public int DecodeBit(int handle)
		{
			int bitID = 1;
			int bitFieldSize = Size();

			while (bitID < bitFieldSize)
			{
				int heapValue = m_buffer[bitID * 2];
				int b = (int)handle < heapValue ? 0 : 1;

				bitID = 2 * bitID | b;
				handle -= heapValue * b;
			}

			return (bitID ^ bitFieldSize);
		}

		/*******************************************************************************
		 * EncodeNode -- Returns the handle associated with the corresponding bitID
		 *
		 * This does the inverse of the DecodeNode routine. Note that this mapping
		 * has the property that any bit set to 0 will be mapped to the ID of the next
		 * bit set to one in the bit field.
		 *
		 */
		public int EncodeBit(int bitID)
		{
			int bitFieldSize = Size();
			int arrayID = bitID + bitFieldSize;
			int handle = 0;

			while (arrayID > 1u)
			{
				uint sibling = (uint)arrayID & (~1u);
				int bitCount = m_buffer[sibling];

				handle += (int)(arrayID & 1u) * bitCount;
				arrayID = arrayID / 2;
			}

			return handle;
		}

		/*******************************************************************************
		 * Reduce -- Sums the 2 elements below the current slot
		 *
		 */
		public void Reduce()
		{
			int depth = (int)FindMSB(Size());

			// iterate over elements atomically
			while (--depth >= 0)
			{
				int minNodeID = 1 << depth;
				int maxNodeID = 2 << depth;

				for (int j = minNodeID; j < maxNodeID; ++j)
				{
					m_buffer[j] = m_buffer[j * 2] + m_buffer[j * 2 + 1];
				}
			}
		}
	}
}



