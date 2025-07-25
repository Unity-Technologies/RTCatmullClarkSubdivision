/* FrustumCulling.glsl - public domain frustum culling GLSL code
    (Created by Jonathan Dupuy 2014.11.13)

*/
//////////////////////////////////////////////////////////////////////////////
//
// Frustum Culling API
//

bool FrustumCullingTest(float4x4 mvpMatrix, float3 bmin, float3 bmax);
bool FrustumCullingTest(in const float4 frustumPlanes[6], float3 bmin, float3 bmax);

//
//
//// end header file /////////////////////////////////////////////////////


// *****************************************************************************
// Frustum Implementation

/**
 * Extract Frustum Planes from MVP Matrix
 *
 * Based on "Fast Extraction of Viewing Frustum Planes from the World-
 * View-Projection Matrix", by Gil Gribb and Klaus Hartmann.
 * This procedure computes the planes of the frustum and normalizes
 * them.
 */
void LoadFrustum(float4x4 mvpMatrix, out float4 planes[6])
{
    for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 2; ++j) {
        planes[i * 2 + j].x = mvpMatrix[0][3] + (j == 0 ? mvpMatrix[0][i] : -mvpMatrix[0][i]);
        planes[i * 2 + j].y = mvpMatrix[1][3] + (j == 0 ? mvpMatrix[1][i] : -mvpMatrix[1][i]);
        planes[i * 2 + j].z = mvpMatrix[2][3] + (j == 0 ? mvpMatrix[2][i] : -mvpMatrix[2][i]);
        planes[i * 2 + j].w = mvpMatrix[3][3] + (j == 0 ? mvpMatrix[3][i] : -mvpMatrix[3][i]);
        planes[i * 2 + j] *= length(planes[i * 2 + j].xyz);
    }
}

/**
 * Negative Vertex of an AABB
 *
 * This procedure computes the negative vertex of an AABB
 * given a normal.
 * See the View Frustum Culling tutorial @ LightHouse3D.com
 * http://www.lighthouse3d.com/tutorials/view-frustum-culling/geometric-approach-testing-boxes-ii/
 */
float3 NegativeVertex(float3 bmin, float3 bmax, float3 n)
{
    bool3 b = bool3(n.x > 0.0, n.y > 0.0, n.z > 0.0);
    return lerp(bmin, bmax, b);
}

/**
 * Frustum-AABB Culling Test
 *
 * This procedure returns true if the AABB is either inside, or in
 * intersection with the frustum, and false otherwise.
 * The test is based on the View Frustum Culling tutorial @ LightHouse3D.com
 * http://www.lighthouse3d.com/tutorials/view-frustum-culling/geometric-approach-testing-boxes-ii/
 */
bool FrustumCullingTest(in const float4 planes[6], float3 bmin, float3 bmax)
{
    float a = 1.0f;

    for (int i = 0; i < 6 && a >= 0.0f; ++i) {
        float3 n = NegativeVertex(bmin, bmax, planes[i].xyz);

        a = dot(float4(n, 1.0f), planes[i]);
    }

    return (a >= 0.0);
}

bool FrustumCullingTest(float4x4 mvpMatrix, float3 bmin, float3 bmax)
{
    float4 planes[6]; LoadFrustum(mvpMatrix, planes);

    return FrustumCullingTest(planes, bmin, bmax);
}




