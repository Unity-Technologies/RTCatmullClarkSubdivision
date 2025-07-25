#ifndef CATMULLCLARKVERTNODES_INCLUDED
#define CATMULLCLARKVERTNODES_INCLUDED

#define CATMULL_CLARK_DISPLACEMENT

#include "../CatmullClark_Gather.cginc"
#include "../ConcurrentBinaryTree/ConcurrentBinaryTree.cginc"
#include "../CatmullClarkTessellation.cginc"


int u_InitializedUniform;
int u_InitializedAdaptive;
float4x4 u_LocalToWorldMatrix;
float4x4 u_WorldToLocalMatrix;
int u_UniformModeDepth;
float4x4 u_ObjectToClipMatrix;
int u_ControlMeshDebugMode;

// Using uniforms here instead of keywords to fit this in a MaterialPropertyBlock and
// avoid asking the user to create multiple material instances for each
int CATMULL_CLARK_ADAPTIVE_MATPROP;
int CATMULL_CLARK_DISPLACEMENT_MATPROP;


// -----------------------------------------------------------------------------
// Catmull Clark Adaptive Vertex Shader
// -----------------------------------------------------------------------------
void CatmullClarkAdaptiveVert_float(float vertexID, out float3 position, out float3 normal, out float3 tangent, out float2 uv)
{
	int intVertexID = (int)round(vertexID) % 3;
	int intInstanceID = (int)round(vertexID) / 3;

	position = float3(0, 0, 0);
	normal = float3(0, 1, 0);
	tangent = float3(1, 0, 0);
	uv = float2(0, 0);

	[branch] if (u_InitializedAdaptive == 1)
	{
		cbt_Node node = cbt_DecodeNode(intInstanceID);
		int subdID = 0;
		cct_Face face = cct_DecodeFace(node, subdID);
		int bisectionDepth = cct_NodeBisectionDepth(node, subdID);
		int ccDepth = 1 + (bisectionDepth >> 1);
		int stride = (ccs_MaxDepth(subdID) - ccDepth) << 1;
		int3 faceHalfedgeIDs = int3(face.halfedgeIDs[0], face.halfedgeIDs[1], face.halfedgeIDs[2]);
		int vertexHalfedgeID = faceHalfedgeIDs[intVertexID];

		float3 vertexNormal = ccs_HalfedgeNormal(subdID, vertexHalfedgeID << stride, ccs_MaxDepth(subdID));
		float3 vertexTangent = ccs_HalfedgeTangent(subdID, vertexHalfedgeID << stride, vertexNormal, ccs_MaxDepth(subdID));
		float2 vertexUV = ccs_HalfedgeVertexUv(subdID, vertexHalfedgeID << stride, ccs_MaxDepth(subdID));
		float3 vertexPosition = ccs_HalfedgeVertexPoint(subdID, vertexHalfedgeID << stride, ccs_MaxDepth(subdID));

		[branch]
		if (CATMULL_CLARK_DISPLACEMENT_MATPROP == 1)
			vertexPosition = ApplyDisplacementMapWithoutCracks(subdID, vertexHalfedgeID, ccDepth, vertexPosition, vertexNormal);

		// Output to ShaderGraph vertex stage
		position = mul(u_LocalToWorldMatrix, float4(vertexPosition, 1.0)).xyz;
		normal = normalize(mul(float4(vertexNormal, 0.0), u_WorldToLocalMatrix).xyz);
		tangent = mul(u_LocalToWorldMatrix, float4(vertexTangent, 0.0)).xyz;
		uv = vertexUV;
	}
}


// -----------------------------------------------------------------------------
// Catmull Clark Uniform Vertex Shader
// -----------------------------------------------------------------------------
void EvenRule(inout int3 x, int b)
{
	if (b == 0)
	{
		x[2] = (x[0] | 2);
		x[1] = (x[0] | 1);
		x[0] = (x[0] | 0);
	}
	else
	{
		x[0] = (x[2] | 2);
		x[1] = (x[2] | 3);
		x[2] = (x[2] | 0);
	}
}

void OddRule(inout int3 x, int b)
{
	if (b == 0)
	{
		x[2] = (x[1] << 2) | 0;
		x[1] = (x[0] << 2) | 2;
		x[0] = (x[0] << 2) | 0;
	}
	else
	{
		x[0] = (x[1] << 2) | 0;
		x[1] = (x[1] << 2) | 2;
		x[2] = (x[2] << 2) | 0;
	}
}

int3 DecodeFace(int meshID, int faceID, int levelOfDetail)
{
	const int baseHalfedgeID = faceID >> levelOfDetail;
	int3 face = int3(
		4 * baseHalfedgeID,
		4 * baseHalfedgeID + 2,
		4 * ccm_HalfedgeNextID(meshID, baseHalfedgeID));
	int pingPong = 0;

	for (int bitID = levelOfDetail - 1; bitID >= 0; --bitID)
	{
		const int bitValue = (faceID >> bitID) & 1;
		if ((pingPong & 1) == 0)
		{
			EvenRule(face, bitValue);
		}
		else
		{
			OddRule(face, bitValue);
		}
		pingPong ^= 1;
	}

	// swap winding for odd faces
	if ((levelOfDetail & 1) == 0)
	{
		const int tmp = face[0];

		face[0] = face[2];
		face[2] = tmp;
	}

	return face;
}

void CatmullClarkUniformVert_float(float vertexID, out float3 position, out float3 normal, out float3 tangent, out float2 uv)
{
	int intVertexID = (int)round(vertexID) % 3;
	int intInstanceID = (int)round(vertexID) / 3;

	position = float3(0, 0, 0);
	normal = float3(0, 1, 0);
	tangent = float3(1, 0, 0);
	uv = float2(0, 0);

	[branch] if (u_InitializedUniform == 1)
	{
		int subdID = 0;
		int ccDepth = 1 + (u_UniformModeDepth >> 1);
		int3 face = DecodeFace(subdID, intInstanceID, u_UniformModeDepth);
		int halfedgeID = face[intVertexID];

		float3 vertexPosition = ccs_HalfedgeVertexPoint(subdID, halfedgeID, ccDepth);
		float3 vertexNormal = ccs_HalfedgeNormal(subdID, halfedgeID, ccDepth);
		float3 vertexTangent = ccs_HalfedgeTangent(subdID, halfedgeID, vertexNormal, ccDepth);
		float2 vertexUV = ccs_HalfedgeVertexUv(subdID, halfedgeID, ccDepth);

		[branch]
		if (CATMULL_CLARK_DISPLACEMENT_MATPROP == 1)
			vertexPosition = ApplyDisplacementMapWithoutCracks(subdID, halfedgeID, ccDepth, vertexPosition, vertexNormal);

		// Output to ShaderGraph vertex stage
		position = mul(u_LocalToWorldMatrix, float4(vertexPosition, 1.0)).xyz;
		normal = normalize(mul(float4(vertexNormal, 0.0), u_WorldToLocalMatrix).xyz);
		tangent = mul(u_LocalToWorldMatrix, float4(vertexTangent, 0.0)).xyz;
		uv = vertexUV;
	}
}


// -----------------------------------------------------------------------------
// Entry Point Uniform/Adaptive
// -----------------------------------------------------------------------------
void CatmullClarkVert_float(float vertexID, out float3 position, out float3 normal, out float3 tangent, out float2 uv)
{
	[branch]
	if (CATMULL_CLARK_ADAPTIVE_MATPROP == 1)
		CatmullClarkAdaptiveVert_float(vertexID, position, normal, tangent, uv);
	else
		CatmullClarkUniformVert_float(vertexID, position, normal, tangent, uv);
}


// -----------------------------------------------------------------------------
// Catmull Clark Debug Vertex Shader
// -----------------------------------------------------------------------------
void CatmullClarkControlMeshVert_float(float vertexID, out float3 position, out float3 normal, out float3 debugColor)
{
	int intVertexID = (int)round(vertexID) % 3;
	int intInstanceID = (int)round(vertexID) / 3;

	position = float3(0, 0, 0);
	normal = float3(0, 1, 0);
	debugColor = float3(1, 0, 0);

	[branch] if (u_InitializedUniform == 1)
	{
		int subdID = 0;
		int cageDepth = 0;
		int ccDepth = 1 + (cageDepth >> 1);

		// Get CCM halfedge positions + face vertex from CCS first level (triangulation of control mesh)
		int halfedgeID = intInstanceID;
		int nextID = ccm_HalfedgeNextID(subdID, halfedgeID);
		int faceID = ccs_VertexToHalfedgeID(subdID, ccm_VertexCount + ccm_HalfedgeFaceID(subdID, halfedgeID), ccDepth);
		float3x3 facePositions = float3x3
		(
			ccm_HalfedgeVertexPoint(subdID, halfedgeID),
			ccm_HalfedgeVertexPoint(subdID, nextID),
			ccs_HalfedgeVertexPoint(subdID, faceID, ccDepth)
		);

		// Output to ShaderGraph vertex stage
		position = mul(u_LocalToWorldMatrix, float4(facePositions[intVertexID], 1.0)).xyz;
		normal = float3(0, 0, 0);

		// DEBUG MODE : WIREFRAME
		if (u_ControlMeshDebugMode == 0)
		{
			// Find wether this face is a triangle or not
			int halfedgeIterator = ccm_HalfedgeNextID(subdID, halfedgeID);
			int halfedgeCount = 1;
			while (halfedgeIterator != halfedgeID && halfedgeCount < 5)
			{
				halfedgeCount++;
				halfedgeIterator = ccm_HalfedgeNextID(subdID, halfedgeIterator);
			}
			bool isTriangle = halfedgeCount == 3;

			// Find wether there is a boundary
			bool isBoundary = ccm_HalfedgeTwinID(subdID, halfedgeID) < 0 ? 1 : 0;

			// Find wether there is a crease
			bool isCrease = ccm_HalfedgeSharpness(subdID, halfedgeID) > 0.0 ? 1.0 : 0.0;

			// Debug colors
			const float3 wireColor = float3(0, 0, 0);
			const float3 triangleColor = float3(0.6, 1.0, 0.75);
			const float3 quadColor = float3(0.96, 0.96, 0.96);
			const float3 boundaryColor = float3(0, 0, 1);
			const float3 creaseColor = float3(1, 0, 0);
			float3 color;

			// The two vertices of the control mesh edge
			if (intVertexID == 0 || intVertexID == 1)
			{
				color = wireColor;
				if (isBoundary == true)
					color = boundaryColor;
				else if (isCrease == true)
					color = creaseColor;
			}
			// The other vertex
			else
			{
				if (isTriangle == true)
					color = triangleColor;
				else
					color = quadColor;
			}
			debugColor = color;
		}
		// DEBUG MODE : UV
		else if (u_ControlMeshDebugMode == 1)
		{
			float3x2 faceUVs = float3x2
			(
				ccm_HalfedgeVertexUv(subdID, halfedgeID),
				ccm_HalfedgeVertexUv(subdID, nextID),
				ccs_HalfedgeVertexUv(subdID, faceID, ccDepth)
			);
			debugColor = float3(faceUVs[intVertexID], 0.0);
		}
	}
}

#endif // CATMULLCLARKVERTNODES_INCLUDED