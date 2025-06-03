struct cc_Halfedge
{
    int twinID;
    int nextID;
    int prevID;
    int faceID;
    int edgeID;
    int vertexID;
    int uvID;
};

struct cc_Halfedge_Quad
{
    int twinID;
    int edgeID;
    int vertexID;
    int uvID;
};

struct cc_Crease
{
    int nextID;
    int prevID;
    float sharpness;
};


// -----------------------------------------------------------------------------
// Buffers
#ifdef CCM_VERTEXTOHALFEDGE_WRITE
RWStructuredBuffer<int> ccm_VertexToHalfedgeBuffer;
#else
StructuredBuffer<int> ccm_VertexToHalfedgeBuffer;
#endif
#ifdef CCM_EDGETOHALFEDGE_WRITE
RWStructuredBuffer<int> ccm_EdgeToHalfedgeBuffer;
#else
StructuredBuffer<int> ccm_EdgeToHalfedgeBuffer;
#endif
#ifdef CCM_FACETOHALFEDGE_WRITE
RWStructuredBuffer<int> ccm_FaceToHalfedgeBuffer;
#else
StructuredBuffer<int> ccm_FaceToHalfedgeBuffer;
#endif
#ifdef CCM_HALFEDGE_WRITE
RWStructuredBuffer<cc_Halfedge> ccm_HalfedgeBuffer;
#else
StructuredBuffer<cc_Halfedge> ccm_HalfedgeBuffer;
#endif
#ifdef CCM_CREASE_WRITE
RWStructuredBuffer<cc_Crease> ccm_CreaseBuffer;
#else
StructuredBuffer<cc_Crease> ccm_CreaseBuffer;
#endif
#ifdef CCM_UV_WRITE
RWStructuredBuffer<float> ccm_UvBuffer;
#else
StructuredBuffer<float> ccm_UvBuffer;
#endif
#ifdef CCM_VERTEX_WRITE
RWStructuredBuffer<float> ccm_VertexPointBuffer;
#else
StructuredBuffer<float> ccm_VertexPointBuffer;
#endif

uniform int ccm_FaceCount;
uniform int ccm_EdgeCount;
uniform int ccm_HalfedgeCount;
uniform int ccm_CreaseCount;
uniform int ccm_VertexCount;
uniform int ccm_UvCount;

#ifdef CCS_VERTEX_WRITE
RWStructuredBuffer<float3> ccs_VertexPointBuffer;
#else
StructuredBuffer<float3> ccs_VertexPointBuffer;
#endif

#ifdef CCS_HALFEDGE_WRITE
RWStructuredBuffer<cc_Halfedge_Quad> ccs_HalfedgeBuffer;
#else
StructuredBuffer<cc_Halfedge_Quad> ccs_HalfedgeBuffer;
#endif

#ifdef CCS_CREASE_WRITE
RWStructuredBuffer<cc_Crease> ccs_CreaseBuffer;
#else
StructuredBuffer<cc_Crease> ccs_CreaseBuffer;
#endif

uniform int ccs_HalfedgeCount;
uniform int ccs_CreaseCount;

// -----------------------------------------------------------------------------
// Displacement
#ifdef CATMULL_CLARK_DISPLACEMENT
Texture2D _Displacement;
float _DisplacementBase;
float _DisplacementOffset;
float _DisplacementAmplitude;
SamplerState linear_repeat_sampler;
#endif

// -----------------------------------------------------------------------------

// counts at a given Catmull-Clark subdivision depth
int ccm_HalfedgeCountAtDepth(const int meshID, int depth);
int ccm_FaceCountAtDepth     (const int meshID, int depth);
int ccm_FaceCountAtDepth_Fast(const int meshID, int depth);
int ccm_EdgeCountAtDepth     (const int meshID, int depth);
int ccm_EdgeCountAtDepth_Fast(const int meshID, int depth);
int ccm_VertexCountAtDepth     (const int meshID, int depth);
int ccm_VertexCountAtDepth_Fast(const int meshID, int depth);

// data-access (O(1))
int ccm_HalfedgeTwinID(const int meshID, int halfedgeID);
int ccm_HalfedgePrevID(const int meshID, int halfedgeID);
int ccm_HalfedgeNextID(const int meshID, int halfedgeID);
int ccm_HalfedgeFaceID(const int meshID, int halfedgeID);
int ccm_HalfedgeEdgeID(const int meshID, int halfedgeID);
int ccm_HalfedgeVertexID(const int meshID, int halfedgeID);
int ccm_HalfedgeUvID(const int meshID, int halfedgeID);
float ccm_HalfedgeSharpnnes(const int meshID, int halfedgeID);
float3 ccm_HalfedgeVertexPoint(const int meshID, int halfedgeID);
float2 ccm_HalfedgeVertexUv(const int meshID, int halfedgeID);
int ccm_CreaseNextID(const int meshID, int edgeID);
int ccm_CreasePrevID(const int meshID, int edgeID);
float ccm_CreaseSharpness(const int meshID, int edgeID);
float3 ccm_VertexPoint(const int meshID, int vertexID);
float2 ccm_Uv(const int meshID, int uvID);
int ccm_HalfedgeNextID_Quad(int halfedgeID);
int ccm_HalfedgePrevID_Quad(int halfedgeID);
int ccm_HalfedgeFaceID_Quad(int halfedgeID);

// (vertex, edge, face) -> halfedge mappings (O(1))
int ccm_VertexToHalfedgeID(const int meshID, int vertexID);
int ccm_EdgeToHalfedgeID(const int meshID, int edgeID);
int ccm_FaceToHalfedgeID(const int meshID, int faceID);
int ccm_FaceToHalfedgeID_Quad(int faceID);

// halfedge remappings (O(1))
int ccm_NextVertexHalfedgeID(const int meshID, int halfedgeID);
int ccm_PrevVertexHalfedgeID(const int meshID, int halfedgeID);

// subd queries
int ccs_MaxDepth(const int subdID);
int ccs_VertexCount(const int subdID);
int ccs_CumulativeFaceCount(const int subdID);
int ccs_CumulativeEdgeCount(const int subdID);
int ccs_CumulativeCreaseCount(const int subdID);
int ccs_CumulativeVertexCount(const int subdID);
int ccs_CumulativeHalfedgeCount(const int subdID);
int ccs_CumulativeFaceCountAtDepth(const int meshID, int depth);
int ccs_CumulativeEdgeCountAtDepth(const int meshID, int depth);
int ccs_CumulativeCreaseCountAtDepth(const int meshID, int depth);
int ccs_CumulativeVertexCountAtDepth(const int meshID, int depth);
int ccs_CumulativeHalfedgeCountAtDepth(const int meshID, int depth);

// O(1) data-access
int ccs_HalfedgeTwinID(const int subdID, int halfedgeID, int depth);
int ccs_HalfedgeNextID(const int subdID, int halfedgeID, int depth);
int ccs_HalfedgePrevID(const int subdID, int halfedgeID, int depth);
int ccs_HalfedgeFaceID(const int meshID, int halfedgeID, int depth);
int ccs_HalfedgeEdgeID(const int meshID, int halfedgeID, int depth);
int ccs_HalfedgeVertexID(const int subdID, int halfedgeID, int depth);
float3 ccs_HalfedgeVertexPoint(const int subdID, int halfedgeID, int depth);
float2 ccs_HalfedgeVertexUv(const int subdID, int halfedgeID, int depth);
int ccs_CreaseNextID_Fast(const int subdID, int edgeID, int depth);
int ccs_CreaseNextID     (const int subdID, int edgeID, int depth);
int ccs_CreasePrevID_Fast(const int subdID, int edgeID, int depth);
int ccs_CreasePrevID     (const int subdID, int edgeID, int depth);
float ccs_CreaseSharpness_Fast(const int subdID, int edgeID, int depth);
float ccs_CreaseSharpness     (const int subdID, int edgeID, int depth);
float ccs_HalfedgeSharpness(const int meshID, int halfedgeID, int depth);
float3 ccs_VertexPoint(const int subdID, int vertexID);

// halfedge remapping
int ccs_NextVertexHalfedgeID(const int subdID, int halfedgeID, int depth);
int ccs_PrevVertexHalfedgeID(const int subdID, int halfedgeID, int depth);

// (vertex, edge) -> halfedge mappings
int ccs_VertexToHalfedgeID(const int subdID, int vertexID, int depth);
int ccs_EdgeToHalfedgeID(const int subdID, int edgeID, int depth);
int ccs_FaceToHalfedgeID(const int subdID, int edgeID, int depth);

// halfedge normal + tangent
float3 ccs_HalfedgeNormal_Fast(const int subdID, int halfedgeIDn, int depth);
float3 ccs_HalfedgeNormal(const int subdID, int halfedgeID, int depth);
void ccs_HalfedgeTangent_Fast(const int subdID, int halfedgeID, int depth, out float3 sdir, out float3 tdir);
float4 ccs_HalfedgeTangent(const int subdID, int halfedgeID, float3 normal, int depth);


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/*******************************************************************************
 * UV Encoding / Decoding routines
 *
 */
float2 cc__DecodeUv(int uvEncoded)
{
    const uint tmp = uint(uvEncoded);
    const float2 uv = float2(
        ((tmp >>  0) & 0xFFFFu) / 65535.0f,
        ((tmp >> 16) & 0xFFFFu) / 65535.0f
    );

    return uv;
}

int cc__EncodeUv(float2 uv)
{
    const uint u = uint(round(uv[0] * 65535.0f));
    const uint v = uint(round(uv[1] * 65535.0f));
    const uint tmp = ((u & 0xFFFFu) | ((v & 0xFFFFu) << 16));

    return int(tmp);
}


/*******************************************************************************
 * FaceCountAtDepth -- Returns the number of faces at a given subdivision depth
 *
 * The number of faces follows the rule
 *          F^{d+1} = H^d
 * Therefore, the number of half edges at a given subdivision depth d>= 0 is
 *          F^d = 4^{d - 1} H^0,
 * where H0 denotes the number of half-edges of the control cage.
 *
 */
int ccm_FaceCountAtDepth_Fast(const int meshID, int depth)
{
    const int H0 = ccm_HalfedgeCount;

    return (H0 << (2 * (depth - 1)));
}

int ccm_FaceCountAtDepth(const int meshID, int depth)
{
    if (depth == 0) {
        return ccm_FaceCount;
    } else {
        return ccm_FaceCountAtDepth_Fast(meshID, depth);
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
int ccm_EdgeCountAtDepth_Fast(const int meshID, int depth)
{
    const int E0 = ccm_EdgeCount;
    const int H0 = ccm_HalfedgeCount;
    const int tmp = ~(0xFFFFFFFF << depth); // (2^d - 1)

    return ((E0 << 1) + (tmp * H0)) << (depth - 1);
}

int ccm_EdgeCountAtDepth(const int meshID, int depth)
{
    if (depth == 0) {
        return ccm_EdgeCount;
    } else {
        return ccm_EdgeCountAtDepth_Fast(meshID, depth);
    }
}


/*******************************************************************************
 * HalfedgeCountAtDepth -- Returns the number of half edges at a given subd depth
 *
 * The number of half edges is multiplied by 4 at each subdivision step.
 * Therefore, the number of half edges at a given subdivision depth d>= 0 is
 *          4^d H0,
 * where H0 denotes the number of half-edges of the control cage.
 *
 */
int ccm_HalfedgeCountAtDepth(const int meshID, int depth)
{
    const int H0 = ccm_HalfedgeCount;

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
int ccm_CreaseCountAtDepth(const int meshID, int depth)
{
    const int C0 = ccm_CreaseCount;

    return C0 << depth;
}


/*******************************************************************************
 * VertexCountAtDepth -- Returns the number of vertices at a given subd depth
 *
 * The number of vertices follows the rule
 *          V^{d+1} = V^d + E^d + F^d
 * For a quad mesh, the number of vertices at a given subdivision depth d>= 0 is
 *          V^d = V0 + (2^{d-1} - 1)E0 + (2^{d-1} - 1)^2F0,
 * where:
 * - V0 denotes the number of vertices of the control cage
 * - E0 denotes the number of edges of the control cage
 * - F0 denotes the number of faces of the control cage
 * Note that since the input mesh may contain non-quad faces, we compute
 * the first subdivision step by hand and then apply the formula.
 *
 */
int ccm_VertexCountAtDepth_Fast(const int meshID, int depth)
{
    const int V0 = ccm_VertexCount;
    const int F0 = ccm_FaceCount;
    const int E0 = ccm_EdgeCount;
    const int H0 = ccm_HalfedgeCount;
    const int F1 = H0;
    const int E1 = 2 * E0 + H0;
    const int V1 = V0 + E0 + F0;
    const int tmp =  ~(0xFFFFFFFF << (depth - 1)); // 2^{d-1} - 1

    return V1 + tmp * (E1 + tmp * F1);
}

int ccm_VertexCountAtDepth(const int meshID, int depth)
{
    if (depth == 0) {
        return ccm_VertexCount;
    } else {
        return ccm_VertexCountAtDepth_Fast(meshID, depth);
    }
}


/*******************************************************************************
 * Halfedge data accessors
 *
 */
cc_Halfedge ccm__Halfedge(const int meshID, int halfedgeID)
{
    return ccm_HalfedgeBuffer[halfedgeID];
}

int ccm_HalfedgeTwinID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).twinID;
}

int ccm_HalfedgeNextID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).nextID;
}

int ccm_HalfedgePrevID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).prevID;
}

int ccm_HalfedgeVertexID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).vertexID;
}

int ccm_HalfedgeUvID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).uvID;
}

int ccm_HalfedgeEdgeID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).edgeID;
}

int ccm_HalfedgeFaceID(const int meshID, int halfedgeID)
{
    return ccm__Halfedge(meshID, halfedgeID).faceID;
}

float ccm_HalfedgeSharpness(const int meshID, int halfedgeID)
{
    return ccm_CreaseSharpness(meshID, ccm_HalfedgeEdgeID(meshID, halfedgeID));
}

float3 ccm_HalfedgeVertexPoint(const int meshID, int halfedgeID)
{
    return ccm_VertexPoint(meshID, ccm_HalfedgeVertexID(meshID, halfedgeID));
}

float2 ccm_HalfedgeVertexUv(const int meshID, int halfedgeID)
{
    return ccm_Uv(meshID, ccm_HalfedgeUvID(meshID, halfedgeID));
}

cc_Crease ccm__Crease(const int meshID, int edgeID)
{
    return ccm_CreaseBuffer[edgeID];
}

int ccm_CreaseNextID(const int meshID, int edgeID)
{
    return ccm__Crease(meshID, edgeID).nextID;
}

int ccm_CreasePrevID(const int meshID, int edgeID)
{
    return ccm__Crease(meshID, edgeID).prevID;
}

float ccm_CreaseSharpness(const int meshID, int edgeID)
{
    return ccm__Crease(meshID, edgeID).sharpness;
}

int ccm_HalfedgeFaceID_Quad(int halfedgeID)
{
    return halfedgeID >> 2;
}


/*******************************************************************************
 * Halfedge Iteration (Quad-only special case)
 *
 */
int ccm__ScrollFaceHalfedgeID_Quad(int halfedgeID, int dir)
{
    const int base = 3;
    const int localID = (halfedgeID & base) + dir;

    return (halfedgeID & ~base) | (localID & base);
}

int ccm_HalfedgeNextID_Quad(int halfedgeID)
{
    return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, +1);
}

int ccm_HalfedgePrevID_Quad(int halfedgeID)
{
    return ccm__ScrollFaceHalfedgeID_Quad(halfedgeID, -1);
}


/*******************************************************************************
 * Vertex queries
 *
 */
float3 ccm_VertexPoint(const int meshID, int vertexID)
{
    const float x = ccm_VertexPointBuffer[3 * vertexID + 0];
    const float y = ccm_VertexPointBuffer[3 * vertexID + 1];
    const float z = ccm_VertexPointBuffer[3 * vertexID + 2];

    return float3(x, y, z);
}

float2 ccm_Uv(const int meshID, int uvID)
{
    const float x = ccm_UvBuffer[2 * uvID + 0];
    const float y = ccm_UvBuffer[2 * uvID + 1];

    return float2(x, y);
}

/*******************************************************************************
 * VertexToHalfedgeID -- Returns a half edge ID that carries a given vertex
 *
 */
int ccm_VertexToHalfedgeID(const int meshID, int vertexID)
{
    return ccm_VertexToHalfedgeBuffer[vertexID];
}


/*******************************************************************************
 * EdgeToHalfedgeID -- Returns a halfedge associated with a given edge
 *
 */
int ccm_EdgeToHalfedgeID(const int meshID, int edgeID)
{
    return ccm_EdgeToHalfedgeBuffer[edgeID];
}


/*******************************************************************************
 * FaceToHalfedgeID -- Returns a halfedge associated with a given face
 *
 */
int ccm_FaceToHalfedgeID(const int meshID, int faceID)
{
    return ccm_FaceToHalfedgeBuffer[faceID];
}

int ccm_FaceToHalfedgeID_Quad(int faceID)
{
    return faceID << 2;
}


/*******************************************************************************
 * Vertex Halfedge Iteration
 *
 */
int ccm_NextVertexHalfedgeID(const int meshID, int halfedgeID)
{
    const int twinID = ccm_HalfedgeTwinID(meshID, halfedgeID);

    return twinID >= 0 ? ccm_HalfedgeNextID(meshID, twinID) : -1;
}

int ccm_PrevVertexHalfedgeID(const int meshID, int halfedgeID)
{
    const int prevID = ccm_HalfedgePrevID(meshID, halfedgeID);

    return ccm_HalfedgeTwinID(meshID, prevID);
}


/*******************************************************************************
 * FaceCountAtDepth -- Returns the accumulated number of faces up to a given subdivision depth
 *
 */
int ccs_CumulativeFaceCountAtDepth(const int meshID, int depth)
{
    return ccs_CumulativeHalfedgeCountAtDepth(meshID, depth) >> 2;
}

int ccs_CumulativeFaceCount(const int subdID)
{
    return ccs_CumulativeFaceCountAtDepth(subdID, ccs_MaxDepth(subdID));
}


/*******************************************************************************
 * EdgeCountAtDepth -- Returns the accumulated number of edges up to a given subdivision depth
 *
 */
int ccs_CumulativeEdgeCountAtDepth(const int meshID, int depth)
{
    const int H0 = ccm_HalfedgeCount;
    const int E0 = ccm_EdgeCount;
    const int H1 = H0 << 2;
    const int E1 = (E0 << 1) + H0;
    const int D = depth;
    const int A = ~(0xFFFFFFFF << D); //  2^{d} - 1

    return (A * (6 * E1 + A * H1 - H1)) / 6;
}

int ccs_CumulativeEdgeCount(const int subdID)
{
    return ccs_CumulativeEdgeCountAtDepth(subdID, ccs_MaxDepth(subdID));
}


/*******************************************************************************
 * HalfedgeCount -- Returns the total number of half edges stored by the subd
 *
 * The number of half edges is multiplied by 4 at each subdivision step.
 * It follows that the number of half-edges is given by the formula
 *    H = H0 x sum_{d=0}^{D} 4^d
 *      = H0 (4^{D+1} - 1) / 3
 * where D denotes the maximum subdivision depth and H0 the number of
 * half edges in the control mesh.
 *
 */
int ccs_CumulativeHalfedgeCountAtDepth(const int meshID, int maxDepth)
{
    const int D = maxDepth;
    const int H0 = ccm_HalfedgeCount;
    const int H1 = H0 << 2;
    const int tmp = ~(0xFFFFFFFF << (D << 1)); // (4^D - 1)

    return (H1 * tmp) / 3;
}

int ccs_CumulativeHalfedgeCount(const int subdID)
{
    return ccs_CumulativeHalfedgeCountAtDepth(subdID, ccs_MaxDepth(subdID));
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
int ccs_CumulativeCreaseCountAtDepth(const int meshID, int maxDepth)
{
    const int D = maxDepth;
    const int C0 = ccm_CreaseCount;
    const int C1 = C0 << 1;
    const int tmp = ~(0xFFFFFFFF << D); // (2^D - 1)

    return (C1 * tmp);
}

int ccs_CumulativeCreaseCount(const int subdID)
{
    return ccs_CumulativeCreaseCountAtDepth(subdID, ccs_MaxDepth(subdID));
}


/*******************************************************************************
 * VertexCount -- Returns the total number of vertices stored by the subd
 *
 * The number of vertices increases according to the following formula at
 * each subdivision step:
 *  Vd+1 = Fd + Ed + Vd
 * It follows that the number of vertices is given by the formula
 *   Vd = d F0 + (2^(d+1) - 1) E0 +
 *      = 4 H0 (4^D - 1) / 3
 * where D denotes the maximum subdivition depth and H0 the number of
 * half edges in the control mesh
 *
 */
int ccs_CumulativeVertexCountAtDepth(const int meshID, int depth)
{
    const int V0 = ccm_VertexCount;
    const int F0 = ccm_FaceCount;
    const int E0 = ccm_EdgeCount;
    const int H0 = ccm_HalfedgeCount;
    const int F1 = H0;
    const int E1 = 2 * E0 + H0;
    const int V1 = V0 + E0 + F0;
    const int D = depth;
    const int A =  ~(0xFFFFFFFF << (D     ));     //  2^{d} - 1
    const int B =  ~(0xFFFFFFFF << (D << 1)) / 3; // (4^{d} - 1) / 3

#if 0 // both formulas are equivalent, not sure which is fastest
    return F1 * (B - (A << 1) + D) + E1 * (A - D) + V1 * D;
#else
    return A * (E1 - (F1 << 1)) + B * F1 + D * (F1 - E1 + V1);
#endif
}

int ccs_CumulativeVertexCount(const int subdID)
{
    return ccs_CumulativeVertexCountAtDepth(subdID, ccs_MaxDepth(subdID));
}


/*******************************************************************************
 * Max Depth Query
 *
 */
int ccs_MaxDepth(const int subdID)
{
    const int Hc = ccs_HalfedgeCount;
    const int H0 = ccm_HalfedgeCount;
    const int H1 = H0 << 2;

    return (firstbithigh((3 * Hc) / H1 + 1) >> 1);
}


/*******************************************************************************
 * VertexCount -- Retrieve the number of vertices
 *
 */
int ccs_VertexCount(const int subdID)
{
    return ccm_VertexCountAtDepth(subdID, ccs_MaxDepth(subdID));
}


/*******************************************************************************
 * Crease queries
 *
 */
cc_Crease ccs__Crease(const int subdID, int edgeID, int depth)
{
    const int stride = ccs_CumulativeCreaseCountAtDepth(subdID, depth - 1);

    return ccs_CreaseBuffer[stride + edgeID];
}

int ccs_CreaseNextID_Fast(const int subdID, int edgeID, int depth)
{
    return ccs__Crease(subdID, edgeID, depth).nextID;
}

int ccs_CreaseNextID(const int subdID, int edgeID, int depth)
{
    const int creaseCount = ccm_CreaseCountAtDepth(subdID, depth);

    if (edgeID < creaseCount) {
        return ccs_CreaseNextID_Fast(subdID, edgeID, depth);
    } else {
        return edgeID;
    }
}

int ccs_CreasePrevID_Fast(const int subdID, int edgeID, int depth)
{
    return ccs__Crease(subdID, edgeID, depth).prevID;
}

int ccs_CreasePrevID(const int subdID, int edgeID, int depth)
{
    const int creaseCount = ccm_CreaseCountAtDepth(subdID, depth);

    if (edgeID < creaseCount) {
        return ccs_CreasePrevID_Fast(subdID, edgeID, depth);
    } else {
        return edgeID;
    }
}

float ccs_CreaseSharpness_Fast(const int subdID, int edgeID, int depth)
{
    return ccs__Crease(subdID, edgeID, depth).sharpness;
}

float ccs_CreaseSharpness(const int subdID, int edgeID, int depth)
{
    const int creaseCount = ccm_CreaseCountAtDepth(subdID, depth);

    if (edgeID < creaseCount) {
        return ccs_CreaseSharpness_Fast(subdID, edgeID, depth);
    } else {
        return 0.0f;
    }
}


/*******************************************************************************
 * Halfedge queries
 *
 */
cc_Halfedge_Quad ccs__Halfedge(const int subdID, int halfedgeID, int depth)
{
    const int stride = ccs_CumulativeHalfedgeCountAtDepth(subdID, depth - 1);

    return ccs_HalfedgeBuffer[stride + halfedgeID];
}

int ccs_HalfedgeTwinID(const int subdID, int halfedgeID, int depth)
{
    return ccs__Halfedge(subdID, halfedgeID, depth).twinID;
}

int ccs_HalfedgeNextID(const int subdID, int halfedgeID, int depth)
{
    return ccm_HalfedgeNextID_Quad(halfedgeID);
}

int ccs_HalfedgePrevID(const int subdID, int halfedgeID, int depth)
{
    return ccm_HalfedgePrevID_Quad(halfedgeID);
}

int ccs_HalfedgeFaceID(const int subdID, int halfedgeID, int depth)
{
    return ccm_HalfedgeFaceID_Quad(halfedgeID);
}

int ccs_HalfedgeEdgeID(const int subdID, int halfedgeID, int depth)
{
    return ccs__Halfedge(subdID, halfedgeID, depth).edgeID;
}

int ccs_HalfedgeVertexID(const int subdID, int halfedgeID, int depth)
{
    return ccs__Halfedge(subdID, halfedgeID, depth).vertexID;
}

float ccs_HalfedgeSharpness(const int subdID, int halfedgeID, int depth)
{
    const int edgeID = ccs_HalfedgeEdgeID(subdID, halfedgeID, depth);

    return ccs_CreaseSharpness(subdID, edgeID, depth);
}

float3 ccs_HalfedgeVertexPoint(const int subdID, int halfedgeID, int depth)
{
    return ccs_VertexPoint(subdID, ccs_HalfedgeVertexID(subdID, halfedgeID, depth));
}

int ccs__HalfedgeUvID(const int subdID, int halfedgeID, int depth)
{
    return ccs__Halfedge(subdID, halfedgeID, depth).uvID;
}

float2 ccs_HalfedgeVertexUv(const int subdID, int halfedgeID, int depth)
{
    return cc__DecodeUv(ccs__HalfedgeUvID(subdID, halfedgeID, depth));
}


/*******************************************************************************
 * Vertex queries
 *
 */
float3 ccs_VertexPoint(const int subdID, int vertexID)
{
	return ccs_VertexPointBuffer[vertexID].xyz;
}


/*******************************************************************************
 * Normal computation
 *
 */
float3 ccs_HalfedgeNormal_Fast(const int subdID, int halfedgeID, int depth)
{
    const int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);
    const int prevID = ccm_HalfedgePrevID_Quad(halfedgeID);
    const float3 v0 = ccs_HalfedgeVertexPoint(subdID, halfedgeID, depth);
    const float3 v1 = ccs_HalfedgeVertexPoint(subdID, prevID    , depth);
    const float3 v2 = ccs_HalfedgeVertexPoint(subdID, nextID    , depth);

    return -normalize(cross(v1 - v0, v2 - v0));
}

float3 ccs_HalfedgeNormal(const int subdID, int halfedgeID, int depth)
{
    const float3 halfedgeNormal = ccs_HalfedgeNormal_Fast(subdID, halfedgeID, depth);
    float3 averageNormal = 0.0.rrr;
    int halfedgeIterator;

	int safeGuard = 0;
    for (halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeID, depth);
         halfedgeIterator >= 0 && halfedgeIterator != halfedgeID && safeGuard < 100;
         halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeIterator, depth))
    {
        averageNormal += ccs_HalfedgeNormal_Fast(subdID, halfedgeIterator, depth);
		safeGuard++;
    }

    if (halfedgeIterator < 0)
        return halfedgeNormal;
    else
        return normalize(halfedgeNormal + averageNormal);
}


/*******************************************************************************
 * Tangent computation
 *
 */
void ccs_HalfedgeTangent_Fast(const int subdID, int halfedgeID, int depth, out float3 sdir, out float3 tdir)
{
	const int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);
	const int prevID = ccm_HalfedgePrevID_Quad(halfedgeID);
	const float3 v0 = ccs_HalfedgeVertexPoint(subdID, halfedgeID, depth);
	const float3 v1 = ccs_HalfedgeVertexPoint(subdID, prevID, depth);
	const float3 v2 = ccs_HalfedgeVertexPoint(subdID, nextID, depth);
	const float2 uv0 = ccs_HalfedgeVertexUv(subdID, halfedgeID, depth);
	const float2 uv1 = ccs_HalfedgeVertexUv(subdID, prevID, depth);
	const float2 uv2 = ccs_HalfedgeVertexUv(subdID, nextID, depth);

	float x0 = v1.x - v0.x;
	float x1 = v2.x - v0.x;
	float y0 = v1.y - v0.y;
	float y1 = v2.y - v0.y;
	float z0 = v1.z - v0.z;
	float z1 = v2.z - v0.z;

	float s0 = uv1.x - uv0.x;
	float s1 = uv2.x - uv0.x;
	float t0 = uv1.y - uv0.y;
	float t1 = uv2.y - uv0.y;

	float r = 1.0 / (s0 * t1 - s1 * t0);
	sdir = float3((t1 * x0 - t0 * x1) * r, (t1 * y0 - t0 * y1) * r, (t1 * z0 - t0 * z1) * r);
	tdir = float3((s0 * x1 - s1 * x0) * r, (s0 * y1 - s1 * y0) * r, (s0 * z1 - s1 * z0) * r);
}

float4 ccs_HalfedgeTangent(const int subdID, int halfedgeID, float3 normal, int depth)
{
	float3 halfedgeTangent1, halfedgeTangent2;
	ccs_HalfedgeTangent_Fast(subdID, halfedgeID, depth, halfedgeTangent1, halfedgeTangent2);
	int uvID = ccs__HalfedgeUvID(subdID, halfedgeID, depth);

	float3 averageTangent1 = 0.0.rrr;
	float3 averageTangent2 = 0.0.rrr;
	int halfedgeIterator;
	int safeGuard = 0;
	for (halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeID, depth);
		halfedgeIterator >= 0 && halfedgeIterator != halfedgeID && safeGuard < 100;
		halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeIterator, depth))
	{
		int uvID2 = ccs__HalfedgeUvID(subdID, halfedgeIterator, depth);
		if (uvID2 == uvID)
		{
			float3 tan1, tan2;
			ccs_HalfedgeTangent_Fast(subdID, halfedgeIterator, depth, tan1, tan2);
			averageTangent1 += tan1;
			averageTangent2 += tan2;
		}
		safeGuard++;
	}

	if (halfedgeIterator < 0)
	{
		float4 tangent = float4(normalize(halfedgeTangent1 - normal * dot(normal, halfedgeTangent1)), // Gram-Schmidt orthogonalize
			(dot(cross(normal, halfedgeTangent1), halfedgeTangent2) < 0.0) ? -1.0 : 1.0); // Handedness
		tangent = min(max(tangent, -1.0), 1.0); // Fix possible NaNs when there was no UV gradient available;
		return tangent;
	}
	else
	{
		averageTangent1 += halfedgeTangent1;
		averageTangent2 += halfedgeTangent2;
		float4 tangent = float4(normalize(averageTangent1 - normal * dot(normal, averageTangent1)), // Gram-Schmidt orthogonalize
			(dot(cross(normal, averageTangent1), averageTangent2) < 0.0) ? -1.0 : 1.0); // Handedness
		tangent = min(max(tangent, -1.0), 1.0); // Fix possible NaNs when there was no UV gradient available;
		return tangent;
	}
}


/*******************************************************************************
 * VertexHalfedge Iteration
 *
 */
int ccs_PrevVertexHalfedgeID(const int subdID, int halfedgeID, int depth)
{
    const int prevID = ccm_HalfedgePrevID_Quad(halfedgeID);

    return ccs_HalfedgeTwinID(subdID, prevID, depth);
}

int ccs_NextVertexHalfedgeID(const int subdID, int halfedgeID, int depth)
{
    const int twinID = ccs_HalfedgeTwinID(subdID, halfedgeID, depth);

    return ccm_HalfedgeNextID_Quad(twinID);
}


/*******************************************************************************
 * Edge to Halfedge Mapping
 *
 * This procedure returns one of the ID of one of the half edge that constitutes
 * the edge. This routine has O(depth) complexity.
 *
 */
int ccs__EdgeToHalfedgeID_First(const int meshID, int edgeID)
{
    const int edgeCount = ccm_EdgeCount;

    if /* [2E, 2E + H) */ (edgeID >= 2 * edgeCount) {
        const int halfedgeID = edgeID - 2 * edgeCount;
        const int nextID = ccm_HalfedgeNextID(meshID, halfedgeID);

        return max(4 * halfedgeID + 1, 4 * nextID + 2);

    } else if /* */ ((edgeID & 1) == 1) {
        const int halfedgeID = ccm_EdgeToHalfedgeID(meshID, edgeID >> 1);
        const int nextID = ccm_HalfedgeNextID(meshID, halfedgeID);

        return 4 * nextID + 3;

    } else /* */ {
        const int halfedgeID = ccm_EdgeToHalfedgeID(meshID, edgeID >> 1);

        return 4 * halfedgeID + 0;
    }
}

int ccs_EdgeToHalfedgeID(const int meshID, int edgeID, int depth)
{
    uint heap = 1u;
    int edgeHalfedgeID = 0;
    int heapDepth = depth;

    // build heap
    for (; heapDepth > 1; --heapDepth) {
        const int edgeCount = ccm_EdgeCountAtDepth_Fast(meshID, heapDepth - 1);

        if /* [2E, 2E + H) */ (edgeID >= 2 * edgeCount) {
            const int halfedgeID = edgeID - 2 * edgeCount;
            const int nextID = ccm_HalfedgeNextID_Quad(halfedgeID);

            edgeHalfedgeID = max(4 * halfedgeID + 1, 4 * nextID + 2);
            break;
        } else {
            heap = (heap << 1) | (edgeID & 1);
            edgeID>>= 1;
        }
    }

    // initialize root cfg
    if (heapDepth == 1) {
        edgeHalfedgeID = ccs__EdgeToHalfedgeID_First(meshID, edgeID);
    }

    // read heap
    while (heap > 1u) {
        if ((heap & 1u) == 1u) {
            const int nextID = ccm_HalfedgeNextID_Quad(edgeHalfedgeID);

            edgeHalfedgeID = 4 * nextID + 3;
        } else {
            edgeHalfedgeID = 4 * edgeHalfedgeID + 0;
        }

        heap>>= 1;
    }

    return edgeHalfedgeID;
}


/*******************************************************************************
 * Vertex to Halfedge Mapping
 *
 * This procedure returns the ID of one of the half edge that connects a
 * given vertex. This routine has O(depth) complexity.
 *
 */
int ccs__VertexToHalfedgeID_First(const int meshID, int vertexID)
{
    const int vertexCount = ccm_VertexCount;
    const int faceCount = ccm_FaceCount;

    if /* [V + F, V + F + E) */ (vertexID >= vertexCount + faceCount) {
        const int edgeID = vertexID - vertexCount - faceCount;

        return 4 * ccm_EdgeToHalfedgeID(meshID, edgeID) + 1;

    } else if /* [V, V + F) */ (vertexID >= vertexCount) {
        const int faceID = vertexID - vertexCount;

        return 4 * ccm_FaceToHalfedgeID(meshID, faceID) + 2;

    } else /* [0, V) */ {

        return 4 * ccm_VertexToHalfedgeID(meshID, vertexID) + 0;
    }
}

int ccs_VertexToHalfedgeID(const int subdID, int vertexID, int depth)
{
    int stride = 0;
    int halfedgeID = 0;
    int heapDepth = depth;

    // build heap
    for (; heapDepth > 1; --heapDepth) {
        const int vertexCount = ccm_VertexCountAtDepth_Fast(subdID, heapDepth - 1);
        const int faceCount = ccm_FaceCountAtDepth_Fast(subdID, heapDepth - 1);

        if /* [V + F, V + F + E) */ (vertexID >= vertexCount + faceCount) {
            const int edgeID = vertexID - faceCount - vertexCount;

            halfedgeID = 4 * ccs_EdgeToHalfedgeID(subdID, edgeID, heapDepth - 1) + 1;
            break;
        } else if /* [V, V + F) */ (vertexID >= vertexCount) {
            const int faceID = vertexID - vertexCount;

            halfedgeID = 4 * ccm_FaceToHalfedgeID_Quad(faceID) + 2;
            break;
        } else /* [0, V) */ {
            stride+= 2;
        }
    }

    // initialize root cfg
    if (heapDepth == 1) {
        halfedgeID = ccs__VertexToHalfedgeID_First(subdID, vertexID);
    }

    return (halfedgeID << stride);
}



#ifdef CCS_VERTEX_WRITE
void WriteVertex(const int cageID, int vertexID, in const float3 vertex)
{
	ccs_VertexPointBuffer[vertexID] = vertex.xyz;
}
#endif



// -----------------------------------------------------------------------------
// Catmull Clark Displacement Mapping
// -----------------------------------------------------------------------------
#ifdef CATMULL_CLARK_DISPLACEMENT
//float3 ApplyDisplacementMapWithoutCracks(int subdID, int vertexHalfedgeID, int depth, in const float3 position, in const float3 normal)
//{
//	int vertexUVID = ccs__HalfedgeUvID(subdID, vertexHalfedgeID, depth);
//
//	// This halfedgeID's displacement sample
//	float2 uv = ccs_HalfedgeVertexUv(subdID, vertexHalfedgeID, depth);
//	float vertexDisplacementSample = _Displacement.SampleLevel(linear_repeat_sampler, uv, 0.0).r;
//
//	// Iterate over neighbouring halfedges to average displacement at UV seams
//	float averageCount = 1.0;
//	float averageDisplacementSample = 0.0;
//	int halfedgeIterator;
//	int safeGuard = 0;
//	for (halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, vertexHalfedgeID, depth);
//		halfedgeIterator >= 0 && halfedgeIterator != vertexHalfedgeID && safeGuard < 100;
//		halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeIterator, depth))
//	{
//		int uvID = ccs__HalfedgeUvID(subdID, halfedgeIterator, depth);
//		float2 otherUV = ccs_HalfedgeVertexUv(subdID, halfedgeIterator, depth);
//		averageDisplacementSample += _Displacement.SampleLevel(linear_repeat_sampler, otherUV, 0.0).r;
//		averageCount += 1.0;
//
//		safeGuard++;
//	}
//
//	float displacement;
//	if (halfedgeIterator < 0)
//		displacement = vertexDisplacementSample;
//	else
//		displacement = (vertexDisplacementSample + averageDisplacementSample) / averageCount;
//
//	displacement = (displacement - _DisplacementBase) * _DisplacementAmplitude + _DisplacementOffset;
//	return position + normal * displacement * 0.01; // HDRP/Lit uses centimeters for these values
//}

float3 ApplyDisplacementMapWithoutCracks(int subdID, int vertexHalfedgeID, int depth, in const float3 position, in const float3 normal)
{
	int vertexUVID = ccs__HalfedgeUvID(subdID, vertexHalfedgeID, depth);

	// This halfedgeID's displacement sample
	float2 uv = ccs_HalfedgeVertexUv(subdID, vertexHalfedgeID, depth);
	float displacement = _Displacement.SampleLevel(linear_repeat_sampler, uv, 0.0).r;

	// Iterate over neighbouring halfedges to average displacement at UV seams
	int minID = vertexUVID;
	int minIDIterator = -1;
	int halfedgeIterator;
	int safeGuard = 0;
	for (halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, vertexHalfedgeID, depth);
		halfedgeIterator >= 0 && halfedgeIterator != vertexHalfedgeID && safeGuard < 100;
		halfedgeIterator = ccs_PrevVertexHalfedgeID(subdID, halfedgeIterator, depth))
	{
		int uvID = ccs__HalfedgeUvID(subdID, halfedgeIterator, depth);
		if (uvID < minID)
		{
			minID = uvID;
			minIDIterator = halfedgeIterator;
		}
		safeGuard++;
	}

	if (minID != vertexUVID)
	{
		float2 otherUV = ccs_HalfedgeVertexUv(subdID, minIDIterator, depth);
		displacement = _Displacement.SampleLevel(linear_repeat_sampler, otherUV, 0.0).r;
	}

	displacement = (displacement - _DisplacementBase) * _DisplacementAmplitude + _DisplacementOffset;
	return position + normal * displacement * 0.01; // HDRP/Lit uses centimeters for these values
}
#endif