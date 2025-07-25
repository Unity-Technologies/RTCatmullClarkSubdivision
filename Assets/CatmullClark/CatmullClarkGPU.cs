using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEngine;
using Unity.Mathematics;
using Unity.Collections;


public static class CatmullClarkGPU
{
	public const int CC_LOCAL_SIZE_X = 256;

	// mesh data-structure
	public class cc_MeshGPU
	{
		public int vertexCount;
		public int uvCount;
		public int halfedgeCount;
		public int edgeCount;
		public int faceCount;
		public ComputeBuffer vertexToHalfedgeIDs;
		public ComputeBuffer edgeToHalfedgeIDs;
		public ComputeBuffer faceToHalfedgeIDs;
		public ComputeBuffer vertexPoints;
		public ComputeBuffer uvs;
		public ComputeBuffer halfedges;
		public ComputeBuffer creases;
	};

	// subd data-structure
	public class cc_SubdGPU
	{
		public int halfedgeCount;
		public int creaseCount;
		public int vertexCount;
		public cc_MeshGPU cage;
		public ComputeBuffer halfedges;
		public ComputeBuffer creases;
		public ComputeBuffer vertexPoints;
		public int maxDepth;
	};

	// Compute shaders
	public static ComputeShader cc_CreasedCageEdgePoints_Gather = (ComputeShader)Resources.Load("cc_CreasedCageEdgePoints_Gather");
	public static ComputeShader cc_CreasedCageFacePoints_Gather = (ComputeShader)Resources.Load("cc_CreasedCageFacePoints_Gather");
	public static ComputeShader cc_CreasedCageVertexPoints_Gather = (ComputeShader)Resources.Load("cc_CreasedCageVertexPoints_Gather");
	public static ComputeShader cc_CreasedEdgePoints_Gather = (ComputeShader)Resources.Load("cc_CreasedEdgePoints_Gather");
	public static ComputeShader cc_CreasedFacePoints_Gather = (ComputeShader)Resources.Load("cc_CreasedFacePoints_Gather");
	public static ComputeShader cc_CreasedVertexPoints_Gather = (ComputeShader)Resources.Load("cc_CreasedVertexPoints_Gather");
	public static ComputeShader cc_RefineCageCreases = (ComputeShader)Resources.Load("cc_RefineCageCreases");
	public static ComputeShader cc_RefineCageHalfedges = (ComputeShader)Resources.Load("cc_RefineCageHalfedges");
	public static ComputeShader cc_RefineCageVertexUvs = (ComputeShader)Resources.Load("cc_RefineCageVertexUvs");
	public static ComputeShader cc_RefineCreases = (ComputeShader)Resources.Load("cc_RefineCreases");
	public static ComputeShader cc_RefineHalfedges = (ComputeShader)Resources.Load("cc_RefineHalfedges");
	public static ComputeShader cc_RefineVertexUvs = (ComputeShader)Resources.Load("cc_RefineVertexUvs");

	public static ComputeShader ccm_BitFieldReduce = (ComputeShader)Resources.Load("ccm_BitFieldReduce");
	public static ComputeShader ccm_LoadFaceMappings = (ComputeShader)Resources.Load("ccm_LoadFaceMappings");
	public static ComputeShader ccm_LoadEdgeMappings = (ComputeShader)Resources.Load("ccm_LoadEdgeMappings");
	public static ComputeShader ccm_ComputeTwins = (ComputeShader)Resources.Load("ccm_ComputeTwins");
	public static ComputeShader ccm_LoadVertexHalfedges = (ComputeShader)Resources.Load("ccm_LoadVertexHalfedges");
	public static ComputeShader ccm_LoadCreases = (ComputeShader)Resources.Load("ccm_LoadCreases");



	/*******************************************************************************
	 * DispatchComputeDistributed -- Dispatch a compute shader by distributing the
	 * thread count to all three dispatch dimensions
	 *
	 */
	private static void DispatchComputeDistributed(ComputeShader compute, int kernel, int threadCount, int groupSizeX)
	{
		int threadGroupCount = (int)math.ceil(threadCount / (float)groupSizeX);
		int threadGroupCountDistributed = (int)math.ceil(math.pow(threadGroupCount, 1.0f / 3.0f));
		compute.SetInt("_DispatchDistributer", threadGroupCountDistributed);
		compute.Dispatch(kernel, threadGroupCountDistributed, threadGroupCountDistributed, threadGroupCountDistributed);
	}




	/*******************************************************************************
	 * FaceCount -- Returns the number of faces
	 *
	 */
	public static int ccm_FaceCount(cc_MeshGPU mesh)
	{
		return mesh.faceCount;
	}


	/*******************************************************************************
	 * EdgeCount -- Returns the number of edges
	 *
	 */
	public static int ccm_EdgeCount(cc_MeshGPU mesh)
	{
		return mesh.edgeCount;
	}


	/*******************************************************************************
	 * CreaseCount -- Returns the number of creases
	 *
	 */
	public static int ccm_CreaseCount(cc_MeshGPU mesh)
	{
		return ccm_EdgeCount(mesh);
	}


	/*******************************************************************************
	 * HalfedgeCount -- Returns the number of halfedges
	 *
	 */
	public static int ccm_HalfedgeCount(cc_MeshGPU mesh)
	{
		return mesh.halfedgeCount;
	}


	/*******************************************************************************
	 * VertexCount -- Returns the number of vertices
	 *
	 */
	public static int ccm_VertexCount(cc_MeshGPU mesh)
	{
		return mesh.vertexCount;
	}


	/*******************************************************************************
	 * UvCount -- Returns the number of uvs
	 *
	 */
	public static int ccm_UvCount(cc_MeshGPU mesh)
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
	public static int ccm_FaceCountAtDepth_Fast(cc_MeshGPU cage, int depth)
	{
		Debug.Assert(depth > 0);
		int H0 = ccm_HalfedgeCount(cage);

		return (H0 << ((depth - 1) << 1));
	}

	public static int ccm_FaceCountAtDepth(cc_MeshGPU cage, int depth)
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
	public static int ccm_EdgeCountAtDepth_Fast(cc_MeshGPU cage, int depth)
	{
		Debug.Assert(depth > 0);
		int E0 = ccm_EdgeCount(cage);
		int H0 = ccm_HalfedgeCount(cage);
		int tmp = (int)(~(0xFFFFFFFF << depth)); // (2^d - 1) // TODO: CHECK UINT->INT CONVERSION

		return ((E0 << 1) + (tmp * H0)) << (depth - 1);
	}

	public static int ccm_EdgeCountAtDepth(cc_MeshGPU cage, int depth)
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
	public static int ccm_HalfedgeCountAtDepth(cc_MeshGPU cage, int depth)
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
	public static int ccm_CreaseCountAtDepth(cc_MeshGPU cage, int depth)
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
	public static int ccm_VertexCountAtDepth_Fast(cc_MeshGPU cage, int depth)
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

	public static int ccm_VertexCountAtDepth(cc_MeshGPU cage, int depth)
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
	public static int ccm_UvCountAtDepth(cc_MeshGPU cage, int depth)
	{
		return ccm_HalfedgeCountAtDepth(cage, depth);
	}


	/*******************************************************************************
	 * Create -- Allocates memory for a mesh of given vertex and halfedge count
	 *
	 */
	public static cc_MeshGPU ccm_Create(
			int vertexCount,
			int uvCount,
			int halfedgeCount,
			int edgeCount,
			int faceCount)
	{
		cc_MeshGPU mesh = new cc_MeshGPU();

		mesh.vertexCount = vertexCount;
		mesh.uvCount = uvCount;
		mesh.halfedgeCount = halfedgeCount;
		mesh.edgeCount = edgeCount;
		mesh.faceCount = faceCount;
		mesh.vertexToHalfedgeIDs = new ComputeBuffer(vertexCount, sizeof(int), ComputeBufferType.Default);
		mesh.edgeToHalfedgeIDs = new ComputeBuffer(edgeCount, sizeof(int), ComputeBufferType.Default);
		mesh.faceToHalfedgeIDs = new ComputeBuffer(faceCount, sizeof(int), ComputeBufferType.Default);
		mesh.halfedges = new ComputeBuffer(halfedgeCount, sizeof(int) * 7, ComputeBufferType.Default);
		mesh.creases = new ComputeBuffer(edgeCount, sizeof(int) * 2 + sizeof(float), ComputeBufferType.Default);
		mesh.vertexPoints = new ComputeBuffer(vertexCount * 3, sizeof(float), ComputeBufferType.Default);
		mesh.uvs = new ComputeBuffer(uvCount > 0 ? uvCount * 2 : 1, sizeof(float), ComputeBufferType.Default);

		return mesh;
	}

	/*******************************************************************************
	 * CreateFromCPU -- Initialize a cc_MeshGPU from a cc_Mesh from CPU library
	 *
	 */
	public static cc_MeshGPU ccm_CreateFromCPU(CatmullClark.cc_Mesh meshCPU)
	{
		cc_MeshGPU meshGPU = ccm_Create(meshCPU.vertexCount,
						  meshCPU.uvCount,
						  meshCPU.halfedgeCount,
						  meshCPU.edgeCount,
						  meshCPU.faceCount);

		meshGPU.vertexToHalfedgeIDs.SetData(meshCPU.vertexToHalfedgeIDs);
		meshGPU.edgeToHalfedgeIDs.SetData(meshCPU.edgeToHalfedgeIDs);
		meshGPU.faceToHalfedgeIDs.SetData(meshCPU.faceToHalfedgeIDs);
		meshGPU.halfedges.SetData(meshCPU.halfedges);
		meshGPU.creases.SetData(meshCPU.creases);
		meshGPU.vertexPoints.SetData(meshCPU.vertexPoints);
		meshGPU.uvs.SetData(meshCPU.uvs);

		return meshGPU;
	}


	/*******************************************************************************
	 * Release -- Releases memory used for a given mesh
	 *
	 */
	public static void ccm_Release(cc_MeshGPU mesh)
	{
		mesh.vertexToHalfedgeIDs.Release();
		mesh.faceToHalfedgeIDs.Release();
		mesh.edgeToHalfedgeIDs.Release();
		mesh.halfedges.Release();
		mesh.creases.Release();
		mesh.vertexPoints.Release();
		mesh.uvs.Release();
	}


	/*******************************************************************************
	 * SetDataToCompute -- Binds mesh data to ComputeShader
	 *
	 */
	public static void ccm_SetDataToCompute(cc_MeshGPU cageGPU, ComputeShader compute, int kernel)
	{
		compute.SetBuffer(kernel, "ccm_VertexToHalfedgeBuffer", cageGPU.vertexToHalfedgeIDs);
		compute.SetBuffer(kernel, "ccm_EdgeToHalfedgeBuffer", cageGPU.edgeToHalfedgeIDs);
		compute.SetBuffer(kernel, "ccm_FaceToHalfedgeBuffer", cageGPU.faceToHalfedgeIDs);
		compute.SetBuffer(kernel, "ccm_HalfedgeBuffer", cageGPU.halfedges);
		compute.SetBuffer(kernel, "ccm_CreaseBuffer", cageGPU.creases);
		compute.SetBuffer(kernel, "ccm_VertexPointBuffer", cageGPU.vertexPoints);
		compute.SetBuffer(kernel, "ccm_UvBuffer", cageGPU.uvs);

		compute.SetInt("ccm_FaceCount", cageGPU.faceCount);
		compute.SetInt("ccm_EdgeCount", cageGPU.edgeCount);
		compute.SetInt("ccm_HalfedgeCount", cageGPU.halfedgeCount);
		compute.SetInt("ccm_CreaseCount", cageGPU.edgeCount);
		compute.SetInt("ccm_VertexCount", cageGPU.vertexCount);
		compute.SetInt("ccm_UvCount", cageGPU.uvCount);
	}

	/*******************************************************************************
	 * SetDataToCompute -- Binds mesh data to Material
	 *
	 */
	public static void ccm_SetDataToMaterial(cc_MeshGPU cageGPU, Material material)
	{
		material.SetBuffer("ccm_VertexToHalfedgeBuffer", cageGPU.vertexToHalfedgeIDs);
		material.SetBuffer("ccm_EdgeToHalfedgeBuffer", cageGPU.edgeToHalfedgeIDs);
		material.SetBuffer("ccm_FaceToHalfedgeBuffer", cageGPU.faceToHalfedgeIDs);
		material.SetBuffer("ccm_HalfedgeBuffer", cageGPU.halfedges);
		material.SetBuffer("ccm_CreaseBuffer", cageGPU.creases);
		material.SetBuffer("ccm_VertexPointBuffer", cageGPU.vertexPoints);
		material.SetBuffer("ccm_UvBuffer", cageGPU.uvs);

		material.SetInt("ccm_FaceCount", cageGPU.faceCount);
		material.SetInt("ccm_EdgeCount", cageGPU.edgeCount);
		material.SetInt("ccm_HalfedgeCount", cageGPU.halfedgeCount);
		material.SetInt("ccm_CreaseCount", cageGPU.edgeCount);
		material.SetInt("ccm_VertexCount", cageGPU.vertexCount);
		material.SetInt("ccm_UvCount", cageGPU.uvCount);
	}

	public static void ccm_SetDataToMaterialPropertyBlock(cc_MeshGPU cageGPU, MaterialPropertyBlock materialPropertyBlock)
	{
		materialPropertyBlock.SetBuffer("ccm_VertexToHalfedgeBuffer", cageGPU.vertexToHalfedgeIDs);
		materialPropertyBlock.SetBuffer("ccm_EdgeToHalfedgeBuffer", cageGPU.edgeToHalfedgeIDs);
		materialPropertyBlock.SetBuffer("ccm_FaceToHalfedgeBuffer", cageGPU.faceToHalfedgeIDs);
		materialPropertyBlock.SetBuffer("ccm_HalfedgeBuffer", cageGPU.halfedges);
		materialPropertyBlock.SetBuffer("ccm_CreaseBuffer", cageGPU.creases);
		materialPropertyBlock.SetBuffer("ccm_VertexPointBuffer", cageGPU.vertexPoints);
		materialPropertyBlock.SetBuffer("ccm_UvBuffer", cageGPU.uvs);

		materialPropertyBlock.SetInt("ccm_FaceCount", cageGPU.faceCount);
		materialPropertyBlock.SetInt("ccm_EdgeCount", cageGPU.edgeCount);
		materialPropertyBlock.SetInt("ccm_HalfedgeCount", cageGPU.halfedgeCount);
		materialPropertyBlock.SetInt("ccm_CreaseCount", cageGPU.edgeCount);
		materialPropertyBlock.SetInt("ccm_VertexCount", cageGPU.vertexCount);
		materialPropertyBlock.SetInt("ccm_UvCount", cageGPU.uvCount);
	}


	/*******************************************************************************
	 * FaceCountAtDepth -- Returns the accumulated number of faces up to a given subdivision depth
	 *
	 */
	public static int ccs_CumulativeFaceCountAtDepth(cc_MeshGPU cage, int depth)
	{
		return ccs_CumulativeHalfedgeCountAtDepth(cage, depth) >> 2;
	}

	public static int ccs_CumulativeFaceCount(cc_SubdGPU subd)
	{
		return ccs_CumulativeFaceCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * EdgeCountAtDepth -- Returns the accumulated number of edges up to a given subdivision depth
	 *
	 */
	public static int ccs_CumulativeEdgeCountAtDepth(cc_MeshGPU cage, int depth)
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

	public static int ccs_CumulativeEdgeCount(cc_SubdGPU subd)
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
	public static int ccs_CumulativeHalfedgeCountAtDepth(cc_MeshGPU cage, int maxDepth)
	{
		Debug.Assert(maxDepth >= 0);
		int D = maxDepth;
		int H0 = ccm_HalfedgeCount(cage);
		int H1 = H0 << 2;
		int tmp = (int)(~(0xFFFFFFFF << (D << 1))); // (4^D - 1) // TODO : CHECK UINT->INT CONVERSION

		return (H1 * tmp) / 3;
	}

	public static int ccs_CumulativeHalfedgeCount(cc_SubdGPU subd)
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
	public static int ccs_CumulativeCreaseCountAtDepth(cc_MeshGPU cage, int maxDepth)
	{
		Debug.Assert(maxDepth >= 0);
		int D = maxDepth;
		int C0 = ccm_CreaseCount(cage);
		int C1 = C0 << 1;
		int tmp = (int)~(0xFFFFFFFF << D); // (2^D - 1) // TODO : CHECK UINT->INT CONVERSION

		return (C1 * tmp);
	}

	public static int ccs_CumulativeCreaseCount(cc_SubdGPU subd)
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
	public static int ccs_CumulativeVertexCountAtDepth(cc_MeshGPU cage, int depth)
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

	public static int ccs_CumulativeVertexCount(cc_SubdGPU subd)
	{
		return ccs_CumulativeVertexCountAtDepth(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * MaxDepth -- Retrieve the maximum subdivision depth of the subd
	 *
	 */
	public static int ccs_MaxDepth(cc_SubdGPU subd)
	{
		return subd.maxDepth;
	}


	/*******************************************************************************
	 * VertexCount -- Retrieve the number of vertices
	 *
	 */
	private static int ccs__VertexCount(cc_MeshGPU cage, int maxDepth)
	{
		return ccm_VertexCountAtDepth_Fast(cage, maxDepth);
	}

	public static int ccs_VertexCount(cc_SubdGPU subd)
	{
		return ccs__VertexCount(subd.cage, ccs_MaxDepth(subd));
	}


	/*******************************************************************************
	 * Create -- Create a subd
	 *
	 */
	public static cc_SubdGPU ccs_Create(int maxDepth, cc_MeshGPU cage)
	{
		cc_SubdGPU subd = new cc_SubdGPU();

		subd.halfedgeCount = ccs_CumulativeHalfedgeCountAtDepth(cage, maxDepth);
		subd.creaseCount = ccs_CumulativeCreaseCountAtDepth(cage, maxDepth);
		subd.vertexCount = ccs__VertexCount(cage, maxDepth);

		subd.maxDepth = maxDepth;
		subd.halfedges = new ComputeBuffer(subd.halfedgeCount, sizeof(int) * 4, ComputeBufferType.Default);
		subd.creases = new ComputeBuffer(subd.creaseCount, sizeof(int) * 2 + sizeof(float), ComputeBufferType.Default);
		subd.vertexPoints = new ComputeBuffer(subd.vertexCount, sizeof(float) * 3, ComputeBufferType.Default);
		subd.cage = cage;

		return subd;
	}


	/*******************************************************************************
	 * CreateFromCPU -- Initialize a cc_SubdGPU from a cc_Subd from CPU library
	 *
	 */
	public static cc_SubdGPU ccs_CreateFromCPU(CatmullClark.cc_Subd subdCPU, cc_MeshGPU meshGPU)
	{
		cc_SubdGPU subdGPU = ccs_Create(subdCPU.maxDepth, meshGPU);

		subdGPU.halfedges.SetData(subdCPU.halfedges);
		subdGPU.creases.SetData(subdCPU.creases);
		subdGPU.vertexPoints.SetData(subdCPU.vertexPoints);

		return subdGPU;
	}


	/*******************************************************************************
	 * Release -- Releases memory used for a given subd
	 *
	 */
	public static void ccs_Release(cc_SubdGPU subd)
	{
		subd.halfedges.Release();
		subd.creases.Release();
		subd.vertexPoints.Release();
	}


	/*******************************************************************************
	 * SetDataToCompute -- Binds subdivision data to ComputeShader
	 *
	 */
	public static void ccs_SetDataToCompute(cc_SubdGPU subdGPU, ComputeShader compute, int kernel)
	{
		compute.SetBuffer(kernel, "ccs_HalfedgeBuffer", subdGPU.halfedges);
		compute.SetBuffer(kernel, "ccs_CreaseBuffer", subdGPU.creases);
		compute.SetBuffer(kernel, "ccs_VertexPointBuffer", subdGPU.vertexPoints);

		compute.SetInt("ccs_HalfedgeCount", subdGPU.halfedgeCount);
		compute.SetInt("ccs_CreaseCount", subdGPU.creaseCount);
	}

	/*******************************************************************************
	 * SetDataToCompute -- Binds subdivision data to Material
	 *
	 */
	public static void ccs_SetDataToMaterial(cc_SubdGPU subdGPU, Material material)
	{
		material.SetBuffer("ccs_HalfedgeBuffer", subdGPU.halfedges);
		material.SetBuffer("ccs_CreaseBuffer", subdGPU.creases);
		material.SetBuffer("ccs_VertexPointBuffer", subdGPU.vertexPoints);

		material.SetInt("ccs_HalfedgeCount", subdGPU.halfedgeCount);
		material.SetInt("ccs_CreaseCount", subdGPU.creaseCount);
	}

	public static void ccs_SetDataToMaterialPropertyBlock(cc_SubdGPU subdGPU, MaterialPropertyBlock materialPropertyBlock)
	{
		materialPropertyBlock.SetBuffer("ccs_HalfedgeBuffer", subdGPU.halfedges);
		materialPropertyBlock.SetBuffer("ccs_CreaseBuffer", subdGPU.creases);
		materialPropertyBlock.SetBuffer("ccs_VertexPointBuffer", subdGPU.vertexPoints);

		materialPropertyBlock.SetInt("ccs_HalfedgeCount", subdGPU.halfedgeCount);
		materialPropertyBlock.SetInt("ccs_CreaseCount", subdGPU.creaseCount);
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
	private static void ccs__CreasedCageFacePoints_Gather(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_CreasedCageFacePoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_FaceCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
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
	private static void ccs__CreasedCageEdgePoints_Gather(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_CreasedCageEdgePoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_EdgeCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
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
	private static void ccs__CreasedCageVertexPoints_Gather(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_CreasedCageVertexPoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_VertexCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
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
	private static void ccs__CreasedFacePoints_Gather(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_CreasedFacePoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_FaceCountAtDepth_Fast(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
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
	private static void ccs__CreasedEdgePoints_Gather(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_CreasedEdgePoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_EdgeCountAtDepth_Fast(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
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
	private static void ccs__CreasedVertexPoints_Gather(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_CreasedVertexPoints_Gather;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_VertexCountAtDepth_Fast(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}



	/*******************************************************************************
	 * RefineVertexPoints -- Computes the result of Catmull Clark subdivision.
	 *
	 */
	public static void ccs_RefineVertexPoints_Gather(cc_SubdGPU subd)
	{
		cc_MeshGPU cage = subd.cage;

		ccs__CreasedCageFacePoints_Gather(subd);
		ccs__CreasedCageEdgePoints_Gather(subd);
		ccs__CreasedCageVertexPoints_Gather(subd);

		for (int depth = 1; depth < ccs_MaxDepth(subd); ++depth)
		{
			ccs__CreasedFacePoints_Gather(subd, depth);
			ccs__CreasedEdgePoints_Gather(subd, depth);
			ccs__CreasedVertexPoints_Gather(subd, depth);
		}
	}


	/*******************************************************************************
	 * RefineCageHalfedges -- Applies halfedge refinement rules on the cage mesh
	 *
	 * This routine computes the halfedges of the control cage after one subdivision
	 * step and stores them in the subd.
	 *
	 */
	private static void ccs__RefineCageHalfedges(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_RefineCageHalfedges;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_HalfedgeCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineHalfedges -- Applies halfedge refinement on the subd
	 *
	 * This routine computes the halfedges of the next subd level.
	 *
	 */
	private static void ccs__RefineHalfedges(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_RefineHalfedges;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_HalfedgeCountAtDepth(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineHalfedges
	 *
	 */
	public static void ccs_RefineHalfedges(cc_SubdGPU subd)
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
	private static void ccs__RefineCageVertexUvs(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_RefineCageVertexUvs;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_HalfedgeCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineVertexUvs -- Applies UV refinement on the subd
	 *
	 * This routine computes the UVs of the next subd level.
	 *
	 */
	private static void ccs__RefineVertexUvs(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_RefineVertexUvs;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_HalfedgeCountAtDepth(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineUvs
	 *
	 */
	public static void ccs_RefineVertexUvs(cc_SubdGPU subd)
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
	private static void ccs__RefineCageCreases(cc_SubdGPU subd)
	{
		ComputeShader compute = cc_RefineCageCreases;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, compute, 0);
		int threadCount = ccm_EdgeCount(subd.cage);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineCreases -- Applies crease subdivision on the subd
	 *
	 * This routine computes the topology of the next subd level.
	 *
	 */
	private static void ccs__RefineCreases(cc_SubdGPU subd, int depth)
	{
		ComputeShader compute = cc_RefineCreases;
		ccm_SetDataToCompute(subd.cage, compute, 0);
		ccs_SetDataToCompute(subd, cc_RefineCreases, 0);
		compute.SetInt("u_Depth", depth);
		int threadCount = ccm_HalfedgeCountAtDepth(subd.cage, depth);
		DispatchComputeDistributed(compute, 0, threadCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * RefineCreases
	 *
	 */
	public static void ccs_RefineCreases(cc_SubdGPU subd)
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
	private static void ccs__RefineTopology(cc_SubdGPU subd)
	{
		ccs_RefineHalfedges(subd);
		ccs_RefineVertexUvs(subd);
		ccs_RefineCreases(subd);
	}

	public static void ccs_Refine_Gather(cc_SubdGPU subd)
	{
		ccs__RefineTopology(subd);
		ccs_RefineVertexPoints_Gather(subd);
	}




	/*******************************************************************************
	 * ComputeTwins -- Computes the twin of each half edge
	 *
	 * This routine is what effectively converts a traditional "indexed mesh"
	 * into a halfedge mesh (in the case where all the primitives are the same).
	 *
	 */
	public static void BitonicSort(cc_MeshGPU mesh, ComputeBuffer hashIDBuffer)
	{
		// Hash ID Buffer needs to have a power of two size (values are padded)
		ccm_SetDataToCompute(mesh, ccm_ComputeTwins, 1);
		ccm_ComputeTwins.SetBuffer(1, "_HashIDBufferRW", hashIDBuffer);

		for (uint d2 = 1; d2 < hashIDBuffer.count; d2 *= 2)
		{
			for (uint d1 = d2; d1 >= 1; d1 /= 2)
			{
				ccm_ComputeTwins.SetInt("_LoopValue1", (int)d1);
				ccm_ComputeTwins.SetInt("_LoopValue2", (int)d2);
				DispatchComputeDistributed(ccm_ComputeTwins, 1, hashIDBuffer.count / 2, 64);
			}
		}
	}

	public static void BitonicSortFast(cc_MeshGPU mesh, ComputeBuffer hashIDBuffer)
	{
		ccm_SetDataToCompute(mesh, ccm_ComputeTwins, 1);
		ccm_SetDataToCompute(mesh, ccm_ComputeTwins, 2);
		ccm_ComputeTwins.SetBuffer(1, "_HashIDBufferRW", hashIDBuffer);
		ccm_ComputeTwins.SetBuffer(2, "_HashIDBufferRW", hashIDBuffer);
		ccm_ComputeTwins.SetInt("_ArraySize", hashIDBuffer.count);

		uint d2;
		for (d2 = 1u; d2 < math.min(128u, hashIDBuffer.count); d2 *= 2u)
		{
			ccm_ComputeTwins.SetInt("_LoopValue1", (int)d2);
			ccm_ComputeTwins.SetInt("_LoopValue2", (int)d2);
			DispatchComputeDistributed(ccm_ComputeTwins, 2, hashIDBuffer.count / 2, 64);
		}

		for (; d2 < hashIDBuffer.count; d2 *= 2u)
		{
			for (uint d1 = d2; d1 >= 128u; d1 /= 2u)
			{
				ccm_ComputeTwins.SetInt("_LoopValue1", (int)d1);
				ccm_ComputeTwins.SetInt("_LoopValue2", (int)d2);
				DispatchComputeDistributed(ccm_ComputeTwins, 1, hashIDBuffer.count / 2, 64);
			}

			ccm_ComputeTwins.SetInt("_LoopValue1", 64);
			ccm_ComputeTwins.SetInt("_LoopValue2", (int)d2);
			DispatchComputeDistributed(ccm_ComputeTwins, 2, hashIDBuffer.count / 2, 64);
		}
	}

	public static void ComputeTwins(cc_MeshGPU mesh)
	{
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		ulong vertexCount = (ulong)ccm_VertexCount(mesh);
		Debug.Assert(vertexCount < 0xFFFFFFFFu, "CatmullClarkMeshRenderer: Unsupported Geometry, too many vertices");

		// Initialize hash ID array
		int nextPowerSize = BitFieldGPU.NextPowerOfTwo(halfedgeCount);
		ComputeBuffer hashIDBuffer = new ComputeBuffer(nextPowerSize, sizeof(uint) * 4);
		ccm_ComputeTwins.SetInt("_ArraySize", hashIDBuffer.count);

		ccm_SetDataToCompute(mesh, ccm_ComputeTwins, 0);
		ccm_ComputeTwins.SetBuffer(0, "_HashIDBufferRW", hashIDBuffer);
		DispatchComputeDistributed(ccm_ComputeTwins, 0, nextPowerSize, CC_LOCAL_SIZE_X);

		// Sort hash IDs with bitonic sort
		BitonicSortFast(mesh, hashIDBuffer);

		// Find potential twin using binary search
		ccm_SetDataToCompute(mesh, ccm_ComputeTwins, 3);
		ccm_ComputeTwins.SetBuffer(3, "_HashIDBuffer", hashIDBuffer);
		DispatchComputeDistributed(ccm_ComputeTwins, 3, halfedgeCount, CC_LOCAL_SIZE_X);

		// Release
		hashIDBuffer.Release();
	}


	/*******************************************************************************
	 * LoadFaceMappings -- Computes the mappings for the faces of the mesh
	 *
	 */
	private static void LoadFaceMappings(cc_MeshGPU mesh, List<MeshTopology> submeshTopologies, List<int[]> submeshIndices)
	{
		// Get submesh topologies and face counts
		int faceCount = 0;
		int[] submeshInfoData = new int[submeshTopologies.Count * 2];
		for (int submeshID = 0; submeshID < submeshTopologies.Count; submeshID++)
		{
			submeshInfoData[submeshID * 2] = submeshTopologies[submeshID] == MeshTopology.Triangles ? 3 : 4;
			submeshInfoData[submeshID * 2 + 1] = submeshIndices[submeshID].Length;
			faceCount += submeshIndices[submeshID].Length / submeshInfoData[submeshID * 2];
		}
		ComputeBuffer submeshInfoBuffer = new ComputeBuffer(submeshInfoData.Length, sizeof(int));
		submeshInfoBuffer.SetData(submeshInfoData);

		mesh.faceCount = faceCount;

		// Initialize bitfield for face mappings
		ccm_SetDataToCompute(mesh, ccm_LoadFaceMappings, 1);
		ccm_LoadFaceMappings.SetInt("_SubmeshCount", submeshTopologies.Count);
		ccm_LoadFaceMappings.SetBuffer(1, "_SubmeshInfoBuffer", submeshInfoBuffer);
		BitFieldGPU faceIterator = new BitFieldGPU(ccm_HalfedgeCount(mesh) + 1);
		faceIterator.BindToCompute(ccm_LoadFaceMappings, 1);
		DispatchComputeDistributed(ccm_LoadFaceMappings, 1, faceCount, CC_LOCAL_SIZE_X);

		// Reduce bitfield
		faceIterator.Reduce();

		// Create halfedge to face mapping
		ccm_SetDataToCompute(mesh, ccm_LoadFaceMappings, 2);
		faceIterator.BindToCompute(ccm_LoadFaceMappings, 2);
		DispatchComputeDistributed(ccm_LoadFaceMappings, 2, ccm_HalfedgeCount(mesh), CC_LOCAL_SIZE_X);

		// Create face to halfedge mapping
		mesh.faceToHalfedgeIDs.Release();
		mesh.faceToHalfedgeIDs = new ComputeBuffer(faceCount, sizeof(int), ComputeBufferType.Default);
		ccm_SetDataToCompute(mesh, ccm_LoadFaceMappings, 2);
		faceIterator.BindToCompute(ccm_LoadFaceMappings, 2);
		DispatchComputeDistributed(ccm_LoadFaceMappings, 2, faceCount, CC_LOCAL_SIZE_X);

		// Release
		submeshInfoBuffer.Release();
	}

	private static void LoadFaceMappingsSimple(cc_MeshGPU mesh, List<MeshTopology> submeshTopologies, List<int[]> submeshIndices)
	{
		// TODO : handle vertex welding for UV discontinuities

		// Get submesh topologies and face counts
		int faceCount = 0;
		int indexCount = 0;
		int[] submeshInfoData = new int[submeshTopologies.Count * 2];
		for (int submeshID = 0; submeshID < submeshTopologies.Count; submeshID++)
		{
			submeshInfoData[submeshID * 2] = submeshTopologies[submeshID] == MeshTopology.Triangles ? 3 : 4;
			submeshInfoData[submeshID * 2 + 1] = submeshIndices[submeshID].Length / submeshInfoData[submeshID * 2];
			faceCount += submeshInfoData[submeshID * 2 + 1];
			indexCount += submeshInfoData[submeshID * 2 + 1] * submeshInfoData[submeshID * 2];
		}
		ComputeBuffer submeshInfoBuffer = new ComputeBuffer(submeshInfoData.Length, sizeof(int));
		submeshInfoBuffer.SetData(submeshInfoData);

		// Build merged index buffer
		NativeArray<int> mergedIndexArray2 = new NativeArray<int>(indexCount, Allocator.Temp);
		int counter = 0;
		for (int submeshID = 0; submeshID < submeshTopologies.Count; submeshID++)
		{
			NativeArray<int>.Copy(submeshIndices[submeshID], 0, mergedIndexArray2, counter, submeshIndices[submeshID].Length);
			counter += submeshIndices[submeshID].Length;
		}
		ComputeBuffer mergedIndexBuffer = new ComputeBuffer(indexCount, sizeof(int));
		mergedIndexBuffer.SetData(mergedIndexArray2);
		mergedIndexArray2.Dispose();

		mesh.faceCount = faceCount;
		mesh.faceToHalfedgeIDs.Release();
		mesh.faceToHalfedgeIDs = new ComputeBuffer(faceCount, sizeof(int), ComputeBufferType.Default);

		// Initialize both mappings at once
		ccm_SetDataToCompute(mesh, ccm_LoadFaceMappings, 0);
		ccm_LoadFaceMappings.SetInt("_SubmeshCount", submeshTopologies.Count);
		ccm_LoadFaceMappings.SetBuffer(0, "_SubmeshInfoBuffer", submeshInfoBuffer);
		ccm_LoadFaceMappings.SetBuffer(0, "_MergedIndexBuffer", mergedIndexBuffer);
		DispatchComputeDistributed(ccm_LoadFaceMappings, 0, faceCount, CC_LOCAL_SIZE_X);

		// Release
		submeshInfoBuffer.Release();
		mergedIndexBuffer.Release();
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
	private static void LoadEdgeMappings(cc_MeshGPU mesh, int precomputedEdgeCount)
	{
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		BitFieldGPU edgeIterator = new BitFieldGPU(halfedgeCount);

		// Init bitfield, one bit per edge
		ccm_SetDataToCompute(mesh, ccm_LoadEdgeMappings, 0);
		edgeIterator.BindToCompute(ccm_LoadEdgeMappings, 0);
		DispatchComputeDistributed(ccm_LoadEdgeMappings, 0, halfedgeCount, CC_LOCAL_SIZE_X);

		// Reduce bitfield
		edgeIterator.Reduce();
		int edgeCount;
		if (precomputedEdgeCount < 0)
			edgeCount = edgeIterator.BitCount();
		else
			edgeCount = precomputedEdgeCount;
		mesh.edgeCount = edgeCount;

		// Create halfedge to edge mapping
		ccm_SetDataToCompute(mesh, ccm_LoadEdgeMappings, 1);
		edgeIterator.BindToCompute(ccm_LoadEdgeMappings, 1);
		DispatchComputeDistributed(ccm_LoadEdgeMappings, 1, halfedgeCount, CC_LOCAL_SIZE_X);

		// Create edge to halfedge mapping
		mesh.edgeToHalfedgeIDs.Release();
		mesh.edgeToHalfedgeIDs = new ComputeBuffer(edgeCount, sizeof(int), ComputeBufferType.Default);
		ccm_SetDataToCompute(mesh, ccm_LoadEdgeMappings, 2);
		edgeIterator.BindToCompute(ccm_LoadEdgeMappings, 2);
		DispatchComputeDistributed(ccm_LoadEdgeMappings, 2, edgeCount, CC_LOCAL_SIZE_X);
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
	private static void LoadVertexHalfedges(cc_MeshGPU mesh, out bool crashed)
	{
		crashed = false;
		int halfedgeCount = ccm_HalfedgeCount(mesh);
		int vertexCount = ccm_VertexCount(mesh);

		// Iterate over halfedges
		ComputeBuffer crashFlagBuffer = new ComputeBuffer(1, sizeof(int));
		//int[] crashFlagData = new int[] { 0 };
		//crashFlagBuffer.SetData(crashFlagData);

		mesh.vertexToHalfedgeIDs.Release();
		mesh.vertexToHalfedgeIDs = new ComputeBuffer(vertexCount, sizeof(int), ComputeBufferType.Default);
		ccm_SetDataToCompute(mesh, ccm_LoadVertexHalfedges, 0);
		ccm_LoadVertexHalfedges.SetBuffer(0, "_CrashFlagRW", crashFlagBuffer);
		DispatchComputeDistributed(ccm_LoadVertexHalfedges, 0, halfedgeCount, CC_LOCAL_SIZE_X);

		// Release
		//crashFlagBuffer.GetData(crashFlagData);
		//crashFlagBuffer.Release();

		// Check for infinite loop
		/*if (crashFlagData[0] > 0)
		{
			crashed = true;
			Debug.LogError("CatmullClarkMeshRenderer: Unsupported geometry, " + crashFlagData[0] + " threads resulted in infinite loops");
			return;
		}*/
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
	private static void InitCreasesAndMakeBoundariesSharp(cc_MeshGPU mesh)
	{
		// Iterate over edges
		int edgeCount = ccm_EdgeCount(mesh);
		mesh.creases.Release();
		mesh.creases = new ComputeBuffer(edgeCount, sizeof(int) * 2 + sizeof(float), ComputeBufferType.Default);
		ccm_SetDataToCompute(mesh, ccm_LoadCreases, 0);
		DispatchComputeDistributed(ccm_LoadCreases, 0, edgeCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * ComputeCreaseNeighbors -- Computes the neighbors of each crease
	 *
	 */
	private static void ComputeCreaseNeighbors(cc_MeshGPU mesh)
	{
		// Iterate over edges
		int edgeCount = ccm_EdgeCount(mesh);
		ccm_SetDataToCompute(mesh, ccm_LoadCreases, 1);
		DispatchComputeDistributed(ccm_LoadCreases, 1, edgeCount, CC_LOCAL_SIZE_X);
	}


	/*******************************************************************************
	 * LoadFromUnityGPU -- Loads a CCM mesh from a Unity mesh using the GPU
	 *
	 * Supports both triangle and quad submeshes (see Keep Quads option). Because
	 * polygons with more than 4 vertices are triangulated by the Unity mesh importer,
	 * the subdivision surface won't produce exactly the same result for these polygons.
	 * Currently all submeshes within the mesh are merged into a single halfedge mesh for
	 * simplicity, which means having a different material per submesh is not supported.
	 * 
	 * UNFINISHED : doesn't support vertex welding, crashes unity
	 */
	public static cc_MeshGPU ccm_LoadFromUnityGPU(Mesh unityMesh, ref int[] unityVertexBufferToCCMWeldedBuffer, bool weldVertices)
	{
		Vector3[] vertices = unityMesh.vertices;
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

		// Init default vertex welding indirection buffer // TODO : add vertex welding on GPU
		unityVertexBufferToCCMWeldedBuffer = new int[vertexCount];
		for (int i = 0; i < vertexCount; i++)
			unityVertexBufferToCCMWeldedBuffer[i] = i;

		// Initialize Catmull Clark Mesh structure, load vertices and texcoords
		cc_MeshGPU meshGPU = ccm_Create(vertexCount, uvCount, halfedgeCount, 1, 1);
		meshGPU.vertexPoints.SetData(vertices);
		meshGPU.uvs.SetData(texcoords0);

		bool crashed = false;
		LoadVertexHalfedges(meshGPU, out crashed);
		if (crashed == true)
			return null;

		InitCreasesAndMakeBoundariesSharp(meshGPU);
		// Load creases here
		ComputeCreaseNeighbors(meshGPU);

		return meshGPU;
	}

	public static void AssertEverythingEqual(cc_MeshGPU mesh0, cc_MeshGPU mesh1)
	{
		// Vertex to Halfedge IDs
		int[] vertexToHalfedgeIDs0 = new int[mesh0.vertexCount];
		mesh0.vertexToHalfedgeIDs.GetData(vertexToHalfedgeIDs0);
		int[] vertexToHalfedgeIDs1 = new int[mesh1.vertexCount];
		mesh1.vertexToHalfedgeIDs.GetData(vertexToHalfedgeIDs1);
		for (int i = 0; i < vertexToHalfedgeIDs0.Length; i++)
		{
			Debug.Assert(vertexToHalfedgeIDs0[i] == vertexToHalfedgeIDs1[i]);
		}

		// Edge to Halfedge IDs
		int[] edgeToHalfedgeIDs0 = new int[mesh0.edgeCount];
		mesh0.edgeToHalfedgeIDs.GetData(edgeToHalfedgeIDs0);
		int[] edgeToHalfedgeIDs1 = new int[mesh1.edgeCount];
		mesh1.edgeToHalfedgeIDs.GetData(edgeToHalfedgeIDs1);
		for (int i = 0; i < edgeToHalfedgeIDs0.Length; i++)
		{
			Debug.Assert(edgeToHalfedgeIDs0[i] == edgeToHalfedgeIDs1[i]);
		}

		// Face to Halfedge IDs
		int[] faceToHalfedgeIDs0 = new int[mesh0.faceCount];
		mesh0.faceToHalfedgeIDs.GetData(faceToHalfedgeIDs0);
		int[] faceToHalfedgeIDs1 = new int[mesh1.faceCount];
		mesh1.faceToHalfedgeIDs.GetData(faceToHalfedgeIDs1);
		for (int i = 0; i < faceToHalfedgeIDs0.Length; i++)
		{
			Debug.Assert(faceToHalfedgeIDs0[i] == faceToHalfedgeIDs1[i]);
		}

		// Halfedges
		CatmullClark.cc_Halfedge[] halfedges0 = new CatmullClark.cc_Halfedge[mesh0.halfedgeCount];
		mesh0.halfedges.GetData(halfedges0);
		CatmullClark.cc_Halfedge[] halfedges1 = new CatmullClark.cc_Halfedge[mesh1.halfedgeCount];
		mesh1.halfedges.GetData(halfedges1);
		for (int i = 0; i < halfedges0.Length; i++)
		{
			if (halfedges0[i].twinID != halfedges1[i].twinID)
			{
				int x = 0;
			}
			Debug.Assert(halfedges0[i].twinID == halfedges1[i].twinID);
			Debug.Assert(halfedges0[i].nextID == halfedges1[i].nextID);
			Debug.Assert(halfedges0[i].prevID == halfedges1[i].prevID);
			Debug.Assert(halfedges0[i].faceID == halfedges1[i].faceID);
			Debug.Assert(halfedges0[i].edgeID == halfedges1[i].edgeID);
			Debug.Assert(halfedges0[i].vertexID == halfedges1[i].vertexID);
			Debug.Assert(halfedges0[i].uvID == halfedges1[i].uvID);
		}

		// Creases
		CatmullClark.cc_Crease[] creases0 = new CatmullClark.cc_Crease[mesh0.edgeCount];
		mesh0.creases.GetData(creases0);
		CatmullClark.cc_Crease[] creases1 = new CatmullClark.cc_Crease[mesh1.edgeCount];
		mesh1.creases.GetData(creases1);
		for (int i = 0; i < creases0.Length; i++)
		{
			/*if (creases0[i].sharpness != creases1[i].sharpness)
			{
				Debug.Log("Crease sharpness at " + i + " : " + creases0[i].sharpness + " != " + creases1[i].sharpness);
			}*/
			Debug.Assert(creases0[i].nextID == creases1[i].nextID);
			Debug.Assert(creases0[i].prevID == creases1[i].prevID);
			Debug.Assert(creases0[i].sharpness == creases1[i].sharpness);
		}

		// Vertex Points
		float3[] vertexPoints0 = new float3[mesh0.vertexCount];
		mesh0.vertexPoints.GetData(vertexPoints0);
		float3[] vertexPoints1 = new float3[mesh1.vertexCount];
		mesh1.vertexPoints.GetData(vertexPoints1);
		for (int i = 0; i < vertexPoints0.Length; i++)
		{
			Debug.Assert(vertexPoints0[i].x == vertexPoints1[i].x);
			Debug.Assert(vertexPoints0[i].y == vertexPoints1[i].y);
			Debug.Assert(vertexPoints0[i].z == vertexPoints1[i].z);
		}

		// UVs
		float2[] uvs0 = new float2[mesh0.uvCount];
		mesh0.uvs.GetData(uvs0);
		float2[] uvs1 = new float2[mesh1.uvCount];
		mesh1.uvs.GetData(uvs1);
		for (int i = 0; i < uvs0.Length; i++)
		{
			/*if (uvs0[i].x == uvs1[i].x || uvs0[i].y == uvs1[i].y)
			{
				Debug.Log("UV at " + i + " : (" + uvs0[i].x + ", " + uvs0[i].y + ") != (" + uvs1[i].x + ", " + uvs1[i].y + ")");
			}*/
			Debug.Assert(uvs0[i].x == uvs1[i].x);
			Debug.Assert(uvs0[i].y == uvs1[i].y);
		}
	}


	/*******************************************************************************
	 * BitFieldGPU class, helper for ccm_LoadFromUnity on the GPU
	 * 
	 */
	private class BitFieldGPU
	{
		private int m_bufferSize;
		private ComputeBuffer m_buffer;


		/*******************************************************************************
		 * Bitfield Constructor
		 *
		 */
		public BitFieldGPU(int size)
		{
			int sizePowerOfTwo = NextPowerOfTwo(size);
			NativeArray<int> temp = new NativeArray<int>(2 * sizePowerOfTwo, Allocator.Temp, NativeArrayOptions.ClearMemory);
			temp[0] = sizePowerOfTwo;

			m_buffer = new ComputeBuffer(2 * sizePowerOfTwo, sizeof(int));
			m_buffer.SetData(temp);
			temp.Dispose();

			m_bufferSize = sizePowerOfTwo;
		}

		public void Release()
		{
			m_buffer.Release();
		}

		/*******************************************************************************
		 * Size -- Returns the size of the bitfield
		 *
		 */
		public int Size()
		{
			return m_bufferSize;
		}

		/*******************************************************************************
		 * NextPowerOfTwo -- Returns the upper power of two value
		 *
		 * if the input is already a power of two, its value is returned.
		 *
		 */
		public static int NextPowerOfTwo(int x)
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
		 * GetBit -- Get a specific bit in the bitfield
		 *
		 */
		public int GetBit(int bitID)
		{
			int[] readback = new int[m_buffer.count];
			m_buffer.GetData(readback);
			int offset = Size();
			/*string txt = "GPU: ";
			for (int i = 0; i < 50; i++)
				txt += readback[offset + i] + ", ";
			Debug.Log(txt);*/
			return readback[offset + bitID];
		}

		/*******************************************************************************
		 * BitCount -- Returns the number of bits set to one in the bit field
		 *
		 */
		public int BitCount()
		{
			int[] readback = new int[m_buffer.count];
			m_buffer.GetData(readback);
			return readback[1];
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
				int threadCount = maxNodeID - minNodeID;

				ccm_BitFieldReduce.SetInt("u_WriteBitStart", minNodeID);
				ccm_BitFieldReduce.SetInt("u_WriteBitCount", threadCount);
				ccm_BitFieldReduce.SetBuffer(0, "u_Bitfield", m_buffer);
				DispatchComputeDistributed(ccm_BitFieldReduce, 0, threadCount, CC_LOCAL_SIZE_X);
			}
		}

		public void BindToCompute(ComputeShader compute, int kernel)
		{
			compute.SetBuffer(kernel, "u_Bitfield", m_buffer);
		}
	}
}
