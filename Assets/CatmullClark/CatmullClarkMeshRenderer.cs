using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.Rendering;
using System.Runtime.InteropServices;
using Unity.Jobs;
using Unity.Collections;
#if UNITY_EDITOR
using UnityEditor;
#endif


[ExecuteAlways]
[RequireComponent(typeof(Renderer))]
public class CatmullClarkMeshRenderer : MonoBehaviour
{
	public enum RenderMode
	{
		AdaptiveTessellation,
		UniformTessellation,
		ControlMesh
	}

	public enum ControlmeshDebugMode
	{
		Wireframe,
		UV
	}

	// Parameters
	[Header("Parameters")]
	public bool resetButton = false;
	public bool pauseAdaptiveTessellation = false;
	public RenderMode renderMode = RenderMode.AdaptiveTessellation;
	public int maxSubdivisionDepth = 2;
	public int adaptiveTargetTrianglePixelLength = 8;
	public bool uniformAutoLOD = false;
	public float uniformAutoLODBias = 1.0f;
	public int uniformAutoLODDepth = 4;
	public int uniformManualDepth = 4;
	public ControlmeshDebugMode debugMode = ControlmeshDebugMode.Wireframe;
	public Camera gameCamera;
	public Material catmullClarkRenderMaterial;
	public Material catmullClarkControlMeshMaterial;
	public Material usedCatmullClarkMaterial;
	public int maxPossibleDepth = 1;
	public bool forceConstantUpdatesInEditMode = false;

	// Displacement parameters
	public bool displacementEnabled = false;
	public Texture2D displacement;
	public float displacementBase;
	public float displacementOffset;
	public float displacementAmplitude;

	// Debug parameters // CUSTOM DEMO
	public bool enableDebugWireframe;
	public bool whiteWireframe;
	public bool useVertexWelding = true;
	public bool loadFromUnityGPU = true;
	public bool bakeCCMFileButton = false;
	public string loadFromCCMFile = "";
	public bool disableSkinnedUpdate = false;

	// Metrics
	public ulong byteSize = 0;
	public string vramUsage = "";
	public int currentTriangleCount = 0;

	// Unity mesh renderer data
	private MeshFilter meshFilter;
	private MeshRenderer meshRenderer;
	private SkinnedMeshRenderer skinnedMeshRenderer;
	private ComputeShader copySkinnedVerticesGPU;

	// Rendering Data
	public Camera targetCamera;
	private bool isRenderingActive = false;
	private RenderParams renderParameters;
	private MaterialPropertyBlock materialParameters;

	// CatmullClark Data
	public CatmullClarkGPU.cc_MeshGPU cageGPU;
	private CatmullClarkGPU.cc_SubdGPU subdivisionGPU;
	private int[] unityVertexBufferToCCMWeldedBuffer;
	private ComputeBuffer unityVertexBufferToCCMWeldedBufferGPU;

	// CatmullClarkTessellation Data
	private ConcurrentBinaryTree cbtTree;
	private ComputeBuffer cbtTreeGPU;
	private ComputeShader cbtSumReduction;
	private ComputeShader cbtDispatcher;
	private ComputeShader cctDispatcher;
	private ComputeShader cctUpdate;
	private ComputeBuffer cbtDispatchArgumentsGPU;
	private GraphicsBuffer cctDrawArgumentsGPU;
	private int splitMergePingPong;

	// Critical variable watch
	private Mesh setMesh;
	private int setMaxDepth = -1;
	public RenderMode setRenderMode = RenderMode.AdaptiveTessellation;
	public Material setMaterial = null;
	private int lastRenderedFrame = -1;
	private int initFromSkinnedWaitFrames = 0;

	// TEMP PERF TEST
	public bool tempTest1 = false;
	public bool tempTest2 = false;
	public bool tempTest3 = false;
	public bool tempTest4 = false;
	private Vector3[] precomputedVertices;
	private List<int[]> precomputedIndices;




	// -----------------------------------------------------------------------------
	// Unity callbacks
	// -----------------------------------------------------------------------------
	void OnEnable()
	{
		meshFilter = GetComponent<MeshFilter>();
		meshRenderer = GetComponent<MeshRenderer>();
		skinnedMeshRenderer = GetComponent<SkinnedMeshRenderer>();

		if (gameCamera == null)
			gameCamera = Camera.main;

		cbtSumReduction = ((ComputeShader)Resources.Load("ConcurrentBinaryTreeSumReduction"));
		cbtDispatcher = ((ComputeShader)Resources.Load("ConcurrentBinaryTreeDispatcher"));
		cctDispatcher = ((ComputeShader)Resources.Load("cct_Dispatcher"));
		cctUpdate = ((ComputeShader)Resources.Load("cct_Update"));

		ChooseSceneOrGameCamera();

		// TEMP PERF TESTMesh sourceMesh;
		Mesh sourceMesh = null;
		if (meshFilter != null)
			sourceMesh = meshFilter.sharedMesh;
		else if (skinnedMeshRenderer != null)
			sourceMesh = skinnedMeshRenderer.sharedMesh;
		precomputedVertices = sourceMesh.vertices;
		precomputedIndices = new List<int[]>();
		for (int submeshID = 0; submeshID < sourceMesh.subMeshCount; submeshID++)
		{
			MeshTopology topology = sourceMesh.GetTopology(submeshID);
			if (topology != MeshTopology.Triangles && topology != MeshTopology.Quads)
				continue;
			precomputedIndices.Add(sourceMesh.GetIndices(submeshID));
		}

		InitEverything();
		CatmullClarkTessellationUpdate();
		EnableRendering();

		if (meshRenderer != null)
			meshRenderer.forceRenderingOff = true;
		if (skinnedMeshRenderer != null)
		{
			skinnedMeshRenderer.forceRenderingOff = true;
			skinnedMeshRenderer.updateWhenOffscreen = true;
		}
	}

	void LateUpdate()
	{
		// If we are in Editor and not currently in Game mode, target SceneViewCamera
		ChooseSceneOrGameCamera();

		Mesh usedMesh;
		Renderer rendererToUse = meshRenderer;
		if (skinnedMeshRenderer != null)
		{
			usedMesh = skinnedMeshRenderer.sharedMesh;
			rendererToUse = skinnedMeshRenderer;
		}
		else
			usedMesh = meshFilter.sharedMesh;

		// If no material is set on this script, disable component
		if (usedCatmullClarkMaterial == null)
		{
			setMaterial = null;
			OnDisable();
			enabled = false;
			return;
		}

		// Reset if needed
		if (resetButton == true || setMesh != usedMesh || setMaxDepth != maxSubdivisionDepth || setMaterial != usedCatmullClarkMaterial || setRenderMode != renderMode)
		{
			bool reloadCCM = setMesh != usedMesh;
			bool reloadCCS = reloadCCM || setMaxDepth != maxSubdivisionDepth;

			DisableRendering();
			ReleaseEverything(reloadCCM, reloadCCS);
			InitEverything(reloadCCM, reloadCCS);
			CatmullClarkTessellationUpdate();
			EnableRendering();

			setMesh = usedMesh;
			setMaxDepth = maxSubdivisionDepth;
			setMaterial = usedCatmullClarkMaterial;
			setRenderMode = renderMode;
			resetButton = false;
		}

		// Update fully Adaptive Tessellation
		if (renderMode == RenderMode.AdaptiveTessellation && pauseAdaptiveTessellation == false)
		{
			Plane[] planes = GeometryUtility.CalculateFrustumPlanes(targetCamera);
			if (GeometryUtility.TestPlanesAABB(planes, rendererToUse.bounds) == true)
			{
				CatmullClarkTessellationUpdate();
			}
		}
		// Update automatic LOD level for Uniform Tessellation
		else if (renderMode == RenderMode.UniformTessellation && uniformAutoLOD == true)
		{
			uniformAutoLODDepth = GetUniformAutoLODDepth(rendererToUse.bounds);
		}

		// Runtime : copy skinned animation vertex buffer from SkinnedMeshRenderer to cc_Mesh.vertexPoints and recompute cc_Subd
		if (skinnedMeshRenderer != null && Application.isPlaying == true && disableSkinnedUpdate == false)
		{
			UpdateVertexPositionsFromSkinnedMeshRenderer();
		}

		// TEMP PERF TESTS
		Mesh sourceMesh;
		if (meshFilter != null)
			sourceMesh = meshFilter.sharedMesh;
		else if (skinnedMeshRenderer != null)
			sourceMesh = skinnedMeshRenderer.sharedMesh;
		else
			return;
		if (tempTest1 == true)
		{
			CatmullClark.cc_Mesh cageCPU;
			cageCPU = CatmullClark.ccm_LoadFromUnity(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding, precomputedVertices, precomputedIndices);
			CatmullClarkGPU.ccm_Release(cageGPU);
			cageGPU = CatmullClarkGPU.ccm_CreateFromCPU(cageCPU);
			CatmullClarkGPU.ccm_SetDataToMaterialPropertyBlock(cageGPU, materialParameters);
			subdivisionGPU.cage = cageGPU;
		}
		else if (tempTest2 == true)
		{
			int precomputedEdgeCount = cageGPU.edgeCount;
			CatmullClarkGPU.ccm_Release(cageGPU);
			cageGPU = CatmullClarkGPU.ccm_LoadFromUnityGPU(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding, precomputedVertices, precomputedIndices, precomputedEdgeCount);
			CatmullClarkGPU.ccm_SetDataToMaterialPropertyBlock(cageGPU, materialParameters);
			subdivisionGPU.cage = cageGPU;
		}
		else if (tempTest3 == true)
		{
			CatmullClarkGPU.ccs_Refine_Gather(subdivisionGPU);
		}

		// Editor : workaround for a correct skinned vertex buffer not being available immediatly following domain reload or application quit,
		// do an init from skinned vertices + CCS subdivision update a couple frames after InitEverything()
		if (skinnedMeshRenderer != null && initFromSkinnedWaitFrames > 0)
		{
			initFromSkinnedWaitFrames--;
			if (initFromSkinnedWaitFrames == 0)
				UpdateVertexPositionsFromSkinnedMeshRenderer();
		}

		// Submit rendering instructions
		if (isRenderingActive == true && tempTest4 == false)
		{
			ScheduleDrawInstructions();
		}

		// Debug : Bake intermediate CCM binary file
		if (bakeCCMFileButton == true)
		{
			bakeCCMFileButton = false;
			if (cageGPU != null)
			{
				string filePath = "Assets/StreamingAssets/BakedCCMs/" + usedMesh.name + ".ccm";
				CatmullClark.ccm_Save(CatmullClark.ccm_CreateFromCPU(cageGPU), filePath);
				loadFromCCMFile = "/StreamingAssets/BakedCCMs/" + usedMesh.name + ".ccm";
			}
		}
	}

	void OnDisable()
	{
		DisableRendering();
		ReleaseEverything();
		if (meshRenderer != null)
			meshRenderer.forceRenderingOff = false;
		if (skinnedMeshRenderer != null)
			skinnedMeshRenderer.forceRenderingOff = false;
	}

#if UNITY_EDITOR
	void OnDrawGizmos()
	{
		// Force continuous MonoBehaviour.Update calls in edit mode
		if (forceConstantUpdatesInEditMode == true && Application.isPlaying == false)
		{
			EditorApplication.QueuePlayerLoopUpdate();
			SceneView.RepaintAll();
		}
	}
#endif



	// -----------------------------------------------------------------------------
	// Management
	// -----------------------------------------------------------------------------
	private void InitEverything(bool initCCM = true, bool initCCS = true)
	{
		if (usedCatmullClarkMaterial == null)
			return;

		// Select source mesh
		Mesh sourceMesh;
		if (meshFilter != null)
			sourceMesh = meshFilter.sharedMesh;
		else if (skinnedMeshRenderer != null)
			sourceMesh = skinnedMeshRenderer.sharedMesh;
		else
			return;

		//System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();
		//sw.Start();

		// Create Catmull Clark Mesh
		if (cageGPU == null || initCCM == true)
		{
			if (cageGPU != null)
				CatmullClarkGPU.ccm_Release(cageGPU);

			if (loadFromUnityGPU == false)
			{
				CatmullClark.cc_Mesh cageCPU;
				if (loadFromCCMFile != "")
					cageCPU = CatmullClark.ccm_Load(Application.dataPath + loadFromCCMFile);
				else
					cageCPU = CatmullClark.ccm_LoadFromUnity(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding);
				if (cageCPU == null)
				{
					ReleaseEverything();
					return;
				}
				cageGPU = CatmullClarkGPU.ccm_CreateFromCPU(cageCPU);
			}
			else
			{
				//cageGPU = CatmullClarkGPU.ccm_LoadFromUnityGPU(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding); // TODO : vertex welding for UV discontinuities
				cageGPU = CatmullClarkGPU.ccm_LoadFromUnityGPU(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding);
			}
		}

		//CatmullClark.cc_Crease[] creases0 = new CatmullClark.cc_Crease[cageGPU.edgeCount];
		//cageGPU.creases.GetData(creases0);
		//sw.Stop();
		//Debug.Log("CCM Load " + (loadFromUnityGPU == true ? "GPU" : "") + " : " + sw.ElapsedMilliseconds + "ms");

		// GPU LOAD PARITY TESTING
		/*{
			CatmullClark.cc_Mesh cage = CatmullClark.ccm_LoadFromUnity(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding);
			CatmullClarkGPU.cc_MeshGPU oldCageGPU = CatmullClarkGPU.ccm_CreateFromCPU(cage);
			CatmullClarkGPU.cc_MeshGPU newCageGPU = CatmullClarkGPU.ccm_LoadFromUnityGPU(sourceMesh, ref unityVertexBufferToCCMWeldedBuffer, useVertexWelding);
			CatmullClarkGPU.AssertEverythingEqual(oldCageGPU, newCageGPU);
			CatmullClarkGPU.ccm_Release(oldCageGPU);
			CatmullClarkGPU.ccm_Release(newCageGPU);
		}*/

		// Create Catmull Clark Subdivision
		if (subdivisionGPU == null || initCCS == true)
		{
			//subdivisionGPU = CatmullClarkGPU.ccs_Create(renderMode == RenderMode.ControlMesh ? 1 : maxSubdivisionDepth, cageGPU);
			subdivisionGPU = CatmullClarkGPU.ccs_Create(maxSubdivisionDepth, cageGPU);
			CatmullClarkGPU.ccs_Refine_Gather(subdivisionGPU);
		}

		// Create Catmull Clark Adaptive Tessellation CBT
		if (renderMode == RenderMode.AdaptiveTessellation)
		{
			cbtTree = CatmullClarkTessellation.cct_Create(subdivisionGPU);
			cbtTreeGPU = new ComputeBuffer(cbtTree.HeapUintSize(), sizeof(uint), ComputeBufferType.Default);
			cbtTreeGPU.SetData(cbtTree.GetHeap());
		}

		cctDrawArgumentsGPU = new GraphicsBuffer(GraphicsBuffer.Target.IndirectArguments, 5, sizeof(int));
		cbtDispatchArgumentsGPU = new ComputeBuffer(4, sizeof(uint), ComputeBufferType.IndirectArguments);

		// Set GPU data to compute shaders
		if (renderMode == RenderMode.AdaptiveTessellation)
		{
			cbtSumReduction = Instantiate(cbtSumReduction);
			cbtDispatcher = Instantiate(cbtDispatcher);
			cctDispatcher = Instantiate(cctDispatcher);
			cctUpdate = Instantiate(cctUpdate);
			SetAllDataToComputeShader(cbtDispatcher, 0);
			cbtDispatcher.SetBuffer(0, "u_CbtDispatchBuffer", cbtDispatchArgumentsGPU);
			SetAllDataToComputeShader(cctUpdate, 0);
			SetAllDataToComputeShader(cctUpdate, 1);
			SetAllDataToComputeShader(cctDispatcher, 0);
			cctDispatcher.SetBuffer(0, "u_CctDispatchBuffer", cctDrawArgumentsGPU);
			SetAllDataToComputeShader(cbtSumReduction, 0);
			SetAllDataToComputeShader(cbtSumReduction, 1);
			SetAllDataToComputeShader(cbtSumReduction, 2);
			UpdateTessellationComputeDisplacementParameters();
		}

		// Set GPU data to render material
		materialParameters = new MaterialPropertyBlock();
		materialParameters.SetInt("CATMULL_CLARK_ADAPTIVE_MATPROP", renderMode == RenderMode.AdaptiveTessellation ? 1 : 0);
		CatmullClarkGPU.ccm_SetDataToMaterialPropertyBlock(cageGPU, materialParameters);
		CatmullClarkGPU.ccs_SetDataToMaterialPropertyBlock(subdivisionGPU, materialParameters);
		materialParameters.SetInt("u_InitializedUniform", 1);
		if (cbtTreeGPU != null)
		{
			materialParameters.SetBuffer("u_CbtBuffer", cbtTreeGPU);
			materialParameters.SetInt("u_InitializedAdaptive", 1);
		}
		renderParameters = new RenderParams(usedCatmullClarkMaterial);
		renderParameters.matProps = materialParameters;
		UpdateMaterialDisplacementParameters();

		// Display VRAM size
		vramUsage = UpdateVideoMemoryUsage();

		// Initialize skinning animation stuff
		if (skinnedMeshRenderer != null)
		{
			copySkinnedVerticesGPU = (ComputeShader)Resources.Load("CatmullClarkRendererCopySkinnedVertices");

			// Prepare GPU indirection map from unity vertex buffer to our welded vertex buffer for skinning copy
			unityVertexBufferToCCMWeldedBufferGPU = new ComputeBuffer(unityVertexBufferToCCMWeldedBuffer.Length, sizeof(int), ComputeBufferType.Default);
			unityVertexBufferToCCMWeldedBufferGPU.SetData(unityVertexBufferToCCMWeldedBuffer);

			// Request compute access for skinned vertex buffer
			skinnedMeshRenderer.vertexBufferTarget |= GraphicsBuffer.Target.Raw;

			// Schedule a skinned vertex buffer copy + CCS subdivision update in a couple frames once Unity has initialized it correctly
			initFromSkinnedWaitFrames = 2;
		}

		// Get maximum possible subdivision depth based on max GPU buffer size possible
		GetMaxPossibleSubdivisionDepth();

		// Force render update in edit mode
		lastRenderedFrame = -1;
	}

	public void ReleaseEverything(bool releaseCCM = true, bool releaseCCS = true)
	{
		if(materialParameters != null)
		{
			materialParameters.SetInt("u_InitializedUniform", 0);
			materialParameters.SetInt("u_InitializedAdaptive", 0);
		}

		if (releaseCCM == true)
		{
			if (cageGPU != null)
				CatmullClarkGPU.ccm_Release(cageGPU);
			cageGPU = null;
		}

		if (releaseCCS == true)
		{
			if (subdivisionGPU != null)
				CatmullClarkGPU.ccs_Release(subdivisionGPU);
			subdivisionGPU = null;
		}

		cbtTree = null;
		if (cbtTreeGPU != null)
			cbtTreeGPU.Release();
		cbtTreeGPU = null;

		if (cctDrawArgumentsGPU != null)
			cctDrawArgumentsGPU.Release();
		cctDrawArgumentsGPU = null;

		UpdateVideoMemoryUsage();
	}

	private void ChooseSceneOrGameCamera()
	{
		targetCamera = gameCamera;
#if UNITY_EDITOR
		if (Application.isPlaying == false
			&& SceneView.lastActiveSceneView != null
			&& SceneView.lastActiveSceneView.camera != null)
		{
			targetCamera = SceneView.lastActiveSceneView.camera;
		}
#endif
	}

	public void GetMaxPossibleSubdivisionDepth()
	{
		if (cageGPU == null)
			maxPossibleDepth = 1;

		int depth = 1;
		int cumulativeHalfedgeCount = CatmullClarkGPU.ccs_CumulativeHalfedgeCountAtDepth(cageGPU, depth) * 4;
		while (cumulativeHalfedgeCount > 0 && cumulativeHalfedgeCount < 536870912)
		{
			depth++;
			cumulativeHalfedgeCount = CatmullClarkGPU.ccs_CumulativeHalfedgeCountAtDepth(cageGPU, depth) * 4;
		}

		maxPossibleDepth = depth - 1;
	}

	private void SetAllDataToComputeShader(ComputeShader compute, int kernel)
	{
		CatmullClarkGPU.ccm_SetDataToCompute(cageGPU, compute, kernel);
		CatmullClarkGPU.ccs_SetDataToCompute(subdivisionGPU, compute, kernel);
		compute.SetBuffer(kernel, "u_CbtBuffer", cbtTreeGPU);
	}

	private string UpdateVideoMemoryUsage()
	{
		if (cageGPU == null || subdivisionGPU == null)
		{
			byteSize = 0;
			return "";
		}

		byteSize =
			  (ulong)cageGPU.vertexToHalfedgeIDs.count * (ulong)cageGPU.vertexToHalfedgeIDs.stride
			+ (ulong)cageGPU.edgeToHalfedgeIDs.count * (ulong)cageGPU.edgeToHalfedgeIDs.stride
			+ (ulong)cageGPU.faceToHalfedgeIDs.count * (ulong)cageGPU.faceToHalfedgeIDs.stride
			+ (ulong)cageGPU.vertexPoints.count * (ulong)cageGPU.vertexPoints.stride
			+ (ulong)cageGPU.uvs.count * (ulong)cageGPU.uvs.stride
			+ (ulong)cageGPU.halfedges.count * (ulong)cageGPU.halfedges.stride
			+ (ulong)cageGPU.creases.count * (ulong)cageGPU.creases.stride
			+ (ulong)subdivisionGPU.halfedges.count * (ulong)subdivisionGPU.halfedges.stride
			+ (ulong)subdivisionGPU.creases.count * (ulong)subdivisionGPU.creases.stride
			+ (ulong)subdivisionGPU.vertexPoints.count * (ulong)subdivisionGPU.vertexPoints.stride;

		if (cbtTreeGPU != null)
			byteSize += (ulong)cbtTreeGPU.count * (ulong)cbtTreeGPU.stride;

		if (false && byteSize >= (1 << 30))
		{
			return string.Format("{0} GiB", byteSize >> 30);
		}
		else if (byteSize >= (1 << 20))
		{
			return string.Format("{0} MiB", byteSize >> 20);
		}
		else if (byteSize >= (1 << 10))
		{
			return string.Format("{0} KiB", byteSize >> 10);
		}
		else
		{
			return string.Format("{0} Bytes", byteSize);
		}
	}

	public void UpdateCurrentTriangleCount()
	{
		if (cageGPU == null || cctDrawArgumentsGPU == null)
		{
			currentTriangleCount = 0;
			return;
		}

		if (renderMode == RenderMode.AdaptiveTessellation)
		{
			AsyncGPUReadback.Request(cctDrawArgumentsGPU, (request) =>
			{
				currentTriangleCount = (request.GetData<int>().ToArray()[0] / 3);
			});
		}
		else if (renderMode == RenderMode.UniformTessellation)
		{
			if (uniformAutoLOD == true)
			{
				currentTriangleCount = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU) << (uniformAutoLODDepth - 1);
			}
			else
			{
				currentTriangleCount = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU) << (uniformManualDepth - 1);
			}
		}
		else
		{
			currentTriangleCount = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU);
		}
	}

	public void UpdateMaterialDisplacementParameters()
	{
		if (displacementEnabled == true)
		{
			materialParameters.SetInt("CATMULL_CLARK_DISPLACEMENT_MATPROP", 1);
			if (displacement != null)
				materialParameters.SetTexture("_Displacement", displacement);
			materialParameters.SetFloat("_DisplacementAmplitude", displacementAmplitude);
			materialParameters.SetFloat("_DisplacementOffset", displacementOffset);
			materialParameters.SetFloat("_DisplacementBase", displacementBase);
		}
		else
		{
			materialParameters.SetInt("CATMULL_CLARK_DISPLACEMENT_MATPROP", 0);
		}
	}

	public void UpdateTessellationComputeDisplacementParameters()
	{
		if (displacementEnabled == true)
		{
			cctUpdate.EnableKeyword("CATMULL_CLARK_DISPLACEMENT");
			if (displacement != null)
			{
				cctUpdate.SetTexture(0, "_Displacement", displacement);
				cctUpdate.SetTexture(1,  "_Displacement", displacement);
			}
			cctUpdate.SetFloat("_DisplacementAmplitude", displacementAmplitude);
			cctUpdate.SetFloat("_DisplacementOffset", displacementOffset);
			cctUpdate.SetFloat("_DisplacementBase", displacementBase);
		}
		else
		{
			cctUpdate.DisableKeyword("CATMULL_CLARK_DISPLACEMENT");
		}
	}



	// -----------------------------------------------------------------------------
	// Rendering
	// -----------------------------------------------------------------------------
	private void ScheduleDrawInstructions()
	{
		// Get active UnityEngine.Renderer
		Renderer rendererToUse = meshRenderer;
		if (skinnedMeshRenderer != null)
			rendererToUse = skinnedMeshRenderer;

		if (cctDrawArgumentsGPU == null || usedCatmullClarkMaterial == null || materialParameters == null || rendererToUse == null || isRenderingActive == false)
			return;

		// If object is animated, we copy the special behaviour of SkinnedMeshRenderer: substitute the object's TRS with the root bone's
		Matrix4x4 localToWorld = transform.localToWorldMatrix;
		Matrix4x4 worldToLocal = transform.worldToLocalMatrix;
		if (skinnedMeshRenderer != null && skinnedMeshRenderer.rootBone != null)
		{
			localToWorld = Matrix4x4.TRS(skinnedMeshRenderer.rootBone.position, skinnedMeshRenderer.rootBone.rotation, Vector3.one);
			worldToLocal = localToWorld.inverse;
		}
		materialParameters.SetMatrix("u_LocalToWorldMatrix", localToWorld);
		materialParameters.SetMatrix("u_WorldToLocalMatrix", worldToLocal);
		materialParameters.SetVector("u_ScreenResolution", new Vector2(targetCamera.pixelWidth, targetCamera.pixelHeight)); // CUSTOM DEMO
		materialParameters.SetFloat("_EnableWireframe", 1.0f); // CUSTOM DEMO

		// Update render parameters
		renderParameters.worldBounds = rendererToUse.bounds;
		renderParameters.motionVectorMode = rendererToUse.motionVectorGenerationMode;
		renderParameters.receiveShadows = rendererToUse.receiveShadows;
		renderParameters.shadowCastingMode = rendererToUse.shadowCastingMode;
		renderParameters.lightProbeUsage = rendererToUse.lightProbeUsage;
		renderParameters.reflectionProbeUsage = rendererToUse.reflectionProbeUsage;
		renderParameters.rendererPriority = rendererToUse.rendererPriority;
		renderParameters.renderingLayerMask = rendererToUse.renderingLayerMask;

		// Update displacement parameters
		UpdateMaterialDisplacementParameters();

		// Rendering
		if (renderMode == RenderMode.AdaptiveTessellation)
		{
			Graphics.RenderPrimitivesIndirect(renderParameters, MeshTopology.Triangles, cctDrawArgumentsGPU);
		}
		else if (renderMode == RenderMode.UniformTessellation)
		{
			int depthToUse = uniformAutoLOD == true ? uniformAutoLODDepth : uniformManualDepth;
			int triCountAtDepth = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU) << (depthToUse - 1);
			int vertexCount = triCountAtDepth * 3;
			materialParameters.SetInt("u_UniformModeDepth", depthToUse - 1);
			Graphics.RenderPrimitives(renderParameters, MeshTopology.Triangles, (int)vertexCount);
		}
		else if (renderMode == RenderMode.UniformTessellation)
		{
			int depthToUse = uniformAutoLOD == true ? uniformAutoLODDepth : uniformManualDepth;
			int triCountAtDepth = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU) << (depthToUse - 1);
			int vertexCount = triCountAtDepth * 3;
			materialParameters.SetInt("u_UniformModeDepth", depthToUse - 1);
			Graphics.RenderPrimitives(renderParameters, MeshTopology.Triangles, (int)vertexCount);
		}
		else // Control Mesh Display
		{
			materialParameters.SetInt("u_ControlMeshDebugMode", (int)debugMode);
			int halfedgeCountCCM = CatmullClarkGPU.ccm_HalfedgeCount(cageGPU);
			int vertexCount = halfedgeCountCCM * 3;
			Graphics.RenderPrimitives(renderParameters, MeshTopology.Triangles, (int)vertexCount);
		}
	}

	public void EnableRendering()
	{
		if (isRenderingActive == false)
		{
			isRenderingActive = true;
		}
	}

	public void DisableRendering()
	{
		if (isRenderingActive == true)
		{
			isRenderingActive = false;
		}
	}

	public int GetUniformAutoLODDepth(Bounds bounds)
	{
		// Frustum culling with current target camera (scene camera when Paused, specified game camera when Playing)
		Plane[] planes = GeometryUtility.CalculateFrustumPlanes(targetCamera);
		if (GeometryUtility.TestPlanesAABB(planes, bounds) == false)
		{
			uniformAutoLODDepth = 1;
		}
		else
		{
			float maxBoundsSize = Mathf.Max(bounds.size.x, Mathf.Max(bounds.size.y, bounds.size.z));
			float cameraNodeDistance = Mathf.Max(0.0f, Vector3.Distance(targetCamera.transform.position, bounds.ClosestPoint(targetCamera.transform.position)));
			float screenSpaceSize = (maxBoundsSize * targetCamera.pixelHeight) / (cameraNodeDistance * 2.0f * Mathf.Tan(targetCamera.fieldOfView * 0.5f * Mathf.Deg2Rad));
			uniformAutoLODDepth = 1 + (int)(Mathf.Clamp01((screenSpaceSize * uniformAutoLODBias) / (float)targetCamera.pixelHeight) * (maxSubdivisionDepth * 2 - 1));
		}
		return uniformAutoLODDepth;
	}



	// -----------------------------------------------------------------------------
	// Adaptive Tessellation with concurrent binary tree
	// -----------------------------------------------------------------------------
	private void CatmullClarkTessellationUpdate()
	{
		if (cbtTree == null || cctDrawArgumentsGPU == null)
			return;

		CbtDispatchPass();
		CbtUpdatePass();
		CbtReductionPass();
		CctDispatchPass();
	}

	private void CbtDispatchPass()
	{
		cbtDispatcher.SetBuffer(0, "u_CbtDispatchBuffer", cbtDispatchArgumentsGPU);
		cbtDispatcher.Dispatch(0, 1, 1, 1);
	}

	private void CbtUpdatePass()
	{
		// If object is animated, we copy the special behaviour of SkinnedMeshRenderer and substitute the object's TRS with the root bone's
		Matrix4x4 localToWorld = transform.localToWorldMatrix;
		if (skinnedMeshRenderer != null && skinnedMeshRenderer.rootBone != null)
			localToWorld = Matrix4x4.TRS(skinnedMeshRenderer.rootBone.position, skinnedMeshRenderer.rootBone.rotation, Vector3.one);

		Matrix4x4 modelMatrix = localToWorld;
		Matrix4x4 viewMatrix = targetCamera.worldToCameraMatrix;
		Matrix4x4 projectionMatrix = GL.GetGPUProjectionMatrix(targetCamera.projectionMatrix, true);
		Matrix4x4 modelViewMatrix = viewMatrix * modelMatrix;

		// Compute LOD factor
		int meshletSubdivisionLevel = 0;
		float tmp = 2.0f * Mathf.Tan(Mathf.Deg2Rad * targetCamera.fieldOfView / 2.0f)
			/ targetCamera.pixelHeight * (1 << meshletSubdivisionLevel)
			* adaptiveTargetTrianglePixelLength;
		float lodFactor = -2.0f * math.log2(tmp) + 2.0f;

		// Update displacement parameters every frame for now
		UpdateTessellationComputeDisplacementParameters();

		// CUSTOM DEMO
		cctUpdate.SetFloat("u_ReverseWindingBackfaceCulling", localToWorld.determinant < 0.0f ? 1.0f : 0.0f);

		cctUpdate.SetMatrix("u_ModelViewMatrix", modelViewMatrix);
		cctUpdate.SetMatrix("u_ModelViewProjectionMatrix", projectionMatrix * modelViewMatrix);
		cctUpdate.SetFloat("u_LodFactor", lodFactor);
		cctUpdate.DispatchIndirect(splitMergePingPong, cbtDispatchArgumentsGPU);

		splitMergePingPong = 1 - splitMergePingPong;
	}

	private void CctDispatchPass()
	{
		cctDispatcher.SetBuffer(0, "u_CctDispatchBuffer", cctDrawArgumentsGPU);
		cctDispatcher.Dispatch(0, 1, 1, 1);
	}

	private void CbtReductionPass()
	{
		int maxDepth = cbtTree.MaxDepth();
		int depth = cbtTree.MaxDepth();
		int kGroupSize = 256;

		if (true)
		{
			// Sum reduction prepass processes 32 nodes in a batch.
			// The Settings component made sure the depth is at least 5.
			int nodeCount = 1 << (depth - 5);
			int numGroups = (nodeCount + kGroupSize - 1) / kGroupSize;

			cbtSumReduction.SetInt("u_PassID", depth);
			cbtSumReduction.Dispatch(0, numGroups, 1, 1);
			depth -= 5;
		}

		while (--depth >= 0)
		{
			int nodeCount = 1 << depth;
			int numGroups = (nodeCount + kGroupSize - 1) / kGroupSize;

			cbtSumReduction.SetInt("u_PassID", depth);

			int writeBitStart = cbtTree.NodeBitID(new ConcurrentBinaryTree.Node(1u << depth, depth));
			int writeBitCountPerGroup = (maxDepth + 1 - depth) * math.min(kGroupSize, nodeCount);
			cbtSumReduction.SetInt("u_WriteBitStart", writeBitStart);
			cbtSumReduction.SetInt("u_WriteBitCountPerGroup", writeBitCountPerGroup);
			bool canCombineWrites = (writeBitStart & 31) == 0 && (writeBitCountPerGroup & 31) == 0;
			int kernel = canCombineWrites ? 1 : 2;
			cbtSumReduction.Dispatch(kernel, numGroups, 1, 1);
		}
	}



	// -----------------------------------------------------------------------------
	// Animation handling
	// -----------------------------------------------------------------------------
	private bool UpdateVertexPositionsFromSkinnedMeshRenderer()
	{
		GraphicsBuffer skinnedVertexBuffer = skinnedMeshRenderer.GetVertexBuffer();
		if (unityVertexBufferToCCMWeldedBuffer == null || skinnedMeshRenderer == null || skinnedVertexBuffer == null)
			return false;

		// Get vertex buffer layout
		int stride = skinnedMeshRenderer.sharedMesh.GetVertexBufferStride(0);
		int offset = skinnedMeshRenderer.sharedMesh.GetVertexAttributeOffset(VertexAttribute.Position);
		int dimension = skinnedMeshRenderer.sharedMesh.GetVertexAttributeDimension(VertexAttribute.Position);

		// Copy unity skinned vertices to our CCM representation on the GPU
		copySkinnedVerticesGPU.SetInt("_VertexBufferStride", stride);
		copySkinnedVerticesGPU.SetInt("_VertexPositionOffset", offset);
		copySkinnedVerticesGPU.SetInt("_VertexPositionDimension", dimension);
		CatmullClarkGPU.ccm_SetDataToCompute(cageGPU, copySkinnedVerticesGPU, 0);
		copySkinnedVerticesGPU.SetBuffer(0, "_UnityVertexBufferToCCMWeldedBuffer", unityVertexBufferToCCMWeldedBufferGPU);
		copySkinnedVerticesGPU.SetBuffer(0, "_UnitySkinnedVertices", skinnedVertexBuffer);
		int threadGroupCount = (int)math.ceil(cageGPU.vertexCount / 256.0f);
		copySkinnedVerticesGPU.Dispatch(0, threadGroupCount, 1, 1);

		// Subdivide new animated vertices
		CatmullClarkGPU.ccs_Refine_Gather(subdivisionGPU);
		return true;
	}
}
