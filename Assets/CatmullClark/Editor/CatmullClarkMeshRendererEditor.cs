using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UnityEngine;
using Unity.Mathematics;
using UnityEngine.Rendering;
#if UNITY_EDITOR
using UnityEditor;
#endif


[CustomEditor(typeof(CatmullClarkMeshRenderer))]
[CanEditMultipleObjects]
[ExecuteAlways]
public class CatmullClarkMeshRendererEditor : Editor
{
	private CatmullClarkMeshRenderer ccmRenderer;

	private SerializedProperty maxSubdivisionDepth;
	private SerializedProperty renderMode;
	private SerializedProperty debugMode;
	private SerializedProperty adaptiveTargetTrianglePixelLength;
	private SerializedProperty uniformAutoLOD;
	private SerializedProperty uniformAutoLODBias;
	private SerializedProperty uniformAutoLODDepth;
	private SerializedProperty uniformManualDepth;
	private SerializedProperty catmullClarkRenderMaterial;
	private SerializedProperty usedCatmullClarkMaterial;
	private SerializedProperty setMaterial;
	private SerializedProperty setRenderMode;
	private SerializedProperty resetButton;
	private SerializedProperty pauseAdaptiveTessellation;
	private SerializedProperty vramUsage;
	private SerializedProperty displacementEnabled;
	private SerializedProperty displacement;
	private SerializedProperty displacementBase;
	private SerializedProperty displacementOffset;
	private SerializedProperty displacementAmplitude;
	private SerializedProperty targetCamera;
	private SerializedProperty gameCamera;
	private SerializedProperty forceConstantUpdatesInEditMode;
	private SerializedProperty enableDebugWireframe;
	private SerializedProperty whiteWireframe;
	private SerializedProperty useVertexWelding;
	private SerializedProperty loadFromUnityGPU;
	private SerializedProperty bakeCCMFileButton;
	private SerializedProperty loadFromCCMFile;



	void OnEnable()
	{
		maxSubdivisionDepth = serializedObject.FindProperty("maxSubdivisionDepth");
		renderMode = serializedObject.FindProperty("renderMode");
		debugMode = serializedObject.FindProperty("debugMode");
		adaptiveTargetTrianglePixelLength = serializedObject.FindProperty("adaptiveTargetTrianglePixelLength");
		uniformAutoLOD = serializedObject.FindProperty("uniformAutoLOD");
		uniformAutoLODBias = serializedObject.FindProperty("uniformAutoLODBias");
		uniformAutoLODDepth = serializedObject.FindProperty("uniformAutoLODDepth");
		uniformManualDepth = serializedObject.FindProperty("uniformManualDepth");
		catmullClarkRenderMaterial = serializedObject.FindProperty("catmullClarkRenderMaterial");
		usedCatmullClarkMaterial = serializedObject.FindProperty("usedCatmullClarkMaterial");
		setMaterial = serializedObject.FindProperty("setMaterial");
		setRenderMode = serializedObject.FindProperty("setRenderMode");
		resetButton = serializedObject.FindProperty("resetButton");
		pauseAdaptiveTessellation = serializedObject.FindProperty("pauseAdaptiveTessellation");
		vramUsage = serializedObject.FindProperty("vramUsage");
		displacementEnabled = serializedObject.FindProperty("displacementEnabled");
		displacement = serializedObject.FindProperty("displacement");
		displacementBase = serializedObject.FindProperty("displacementBase");
		displacementOffset = serializedObject.FindProperty("displacementOffset");
		displacementAmplitude = serializedObject.FindProperty("displacementAmplitude");
		targetCamera = serializedObject.FindProperty("targetCamera");
		gameCamera = serializedObject.FindProperty("gameCamera");
		forceConstantUpdatesInEditMode = serializedObject.FindProperty("forceConstantUpdatesInEditMode");
		enableDebugWireframe = serializedObject.FindProperty("enableDebugWireframe");
		whiteWireframe = serializedObject.FindProperty("whiteWireframe");
		useVertexWelding = serializedObject.FindProperty("useVertexWelding");
		loadFromUnityGPU = serializedObject.FindProperty("loadFromUnityGPU");
		bakeCCMFileButton = serializedObject.FindProperty("bakeCCMFileButton");
		loadFromCCMFile = serializedObject.FindProperty("loadFromCCMFile");
	}

	public override bool RequiresConstantRepaint()
	{
		CatmullClarkMeshRenderer[] ccmRenderers = new CatmullClarkMeshRenderer[targets.Length];
		for (int i = 0; i < targets.Length; i++)
			ccmRenderers[i] = (CatmullClarkMeshRenderer)(targets[i]);

		for (int i = 0; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.AdaptiveTessellation && ccmRenderers[i].pauseAdaptiveTessellation == false)
				return true;
		return false;
	}

	public override void OnInspectorGUI()
	{
		serializedObject.Update();

		ccmRenderer = (CatmullClarkMeshRenderer)target;
		CatmullClarkMeshRenderer[] ccmRenderers = new CatmullClarkMeshRenderer[targets.Length];
		for (int i = 0; i < targets.Length; i++)
			ccmRenderers[i] = (CatmullClarkMeshRenderer)(targets[i]);



		// Common Parameters
		EditorGUILayout.Space();
		EditorGUILayout.LabelField("Common Parameters", EditorStyles.boldLabel);

		EditorGUILayout.PropertyField(forceConstantUpdatesInEditMode, new GUIContent("Force Editor Constant Updates"));
		EditorGUILayout.PropertyField(catmullClarkRenderMaterial, new GUIContent("Catmull Clark Material"));
		if (MultiEditingAllSetRenderMaterials(ccmRenderers) == false)
		{
			if (ccmRenderers.Length > 1)
				EditorGUILayout.HelpBox("One or more CatmullClarkMeshRnderers don't have a material set, assign a CatmullClarkVert ShaderGraph material", MessageType.Warning);
			else
				EditorGUILayout.HelpBox("CatmullClarkVert ShaderGraph needed, assign a compatible material then re-enable this component", MessageType.Warning);
		}
		for (int i = 0; i < ccmRenderers.Length; i++)
			ccmRenderers[i].usedCatmullClarkMaterial = ccmRenderers[i].catmullClarkRenderMaterial;

		int minimumMaxDepth = MultiEditingMinimumMaxDepth(ccmRenderers);
		PropertyIntSlider(maxSubdivisionDepth, 1, minimumMaxDepth, new GUIContent("Subdivisions"));



		// Performance metrics
		for (int i = 0; i < ccmRenderers.Length; i++)
			ccmRenderers[i].UpdateCurrentTriangleCount();
		string perfTextBox = "<b>";
		if (ccmRenderers.Length > 1)
			perfTextBox += "Selection ";
		perfTextBox += "VRAM Usage:</b>";
		perfTextBox += "\t\t<i>" + MultiEditingTotalMemoryConsumption(ccmRenderers) + "</i>";
		perfTextBox += "\n<b>";
		if (ccmRenderers.Length > 1)
			perfTextBox += "Selection ";
		perfTextBox += "Triangle Count:</b>";
		perfTextBox += "\t\t<i>" + string.Format("{0:n0}", MultiEditingTotalTriangleCount(ccmRenderers)) + "</i>";
		GUIStyle style1 = new GUIStyle(EditorStyles.centeredGreyMiniLabel);
		style1.alignment = TextAnchor.MiddleLeft;
		EditorGUILayout.LabelField("Performance metrics", style1);
		GUIStyle style2 = new GUIStyle(EditorStyles.helpBox);
		style2.richText = true;
		EditorGUILayout.TextArea(perfTextBox, style2);



		// Rendering Parameters
		EditorGUILayout.Space();
		EditorGUILayout.Space();
		EditorGUILayout.LabelField("Rendering Parameters", EditorStyles.boldLabel);

		EditorGUILayout.PropertyField(renderMode, new GUIContent("Rendering Mode"));

		// Display RenderMode specific parameters only when render mode is identical across selection
		if (MultiEditingAllSameRenderMode(ccmRenderers) == true)
		{
			CatmullClarkMeshRenderer.RenderMode selectionRenderMode = ccmRenderers[0].renderMode;
			if (selectionRenderMode == CatmullClarkMeshRenderer.RenderMode.AdaptiveTessellation)
			{
				EditorGUILayout.PropertyField(pauseAdaptiveTessellation, new GUIContent("Pause Adaptive Tessellation"));
				PropertyIntSlider(adaptiveTargetTrianglePixelLength, 1, 100, new GUIContent("Target Triangle Pixel Length"));
			}
			else if (selectionRenderMode == CatmullClarkMeshRenderer.RenderMode.UniformTessellation)
			{
				EditorGUILayout.BeginHorizontal();
				EditorGUILayout.PropertyField(uniformAutoLOD, new GUIContent("Automatic LOD Selection"));
				if (ccmRenderers.Length == 1 && ccmRenderer.uniformAutoLOD == true)
					EditorGUILayout.LabelField("Current LOD: " + ccmRenderer.uniformAutoLODDepth);
				EditorGUILayout.EndHorizontal();

				// Display auto/manual uniform LOD specific parameters only when autoLOD mode is identical across selection
				if (MultiEditingAllSameAutoLODMode(ccmRenderers) == true)
				{
					bool uniformAutoLOD = ccmRenderers[0].uniformAutoLOD;
					int minSubdivDepth = MultiEditingMinimumSubdivisionDepth(ccmRenderers);
					if (uniformAutoLOD == false)
						PropertyIntSlider(uniformManualDepth, 1, minSubdivDepth * 2, new GUIContent("Level of Detail"));
					else
						PropertyFloatSlider(uniformAutoLODBias, 0.0f, 2.0f, new GUIContent("LOD Bias"));
				}
			}
			else
			{
				// Init control mesh material internally
				for (int i = 0; i < ccmRenderers.Length; i++)
				{
					if (ccmRenderers[i].catmullClarkControlMeshMaterial == null)
					{
						//ccmRenderers[i].catmullClarkControlMeshMaterial = new Material(Shader.Find("Shader Graphs/CatmullClarkControlMeshDebug"));
						//ccmRenderers[i].catmullClarkControlMeshMaterial = new Material(Shader.Find("Custom/URPLitCatmullClarkWireframeControl")); // CUSTOM DEMO
					}
					//ccmRenderers[i].usedCatmullClarkMaterial = ccmRenderers[i].catmullClarkControlMeshMaterial;
				}

				EditorGUILayout.PropertyField(debugMode, new GUIContent("Control Mesh Debug Mode"));
				if (MultiEditingOneWireframeMode(ccmRenderers) == true)
				{
					EditorGUILayout.LabelField("Colour legend", style1);
					EditorGUILayout.TextArea("White polygons: quads\nGreen polygons: triangles\nBlue edges: boundaries\nRed edges: creases", style2);
				}
			}

			// Target camera
			if (MultiEditingAllAdaptive(ccmRenderers) == true)
			{
				EditorGUILayout.PropertyField(gameCamera, new GUIContent("Play Mode Target Camera"));
				GUI.enabled = false;
				EditorGUILayout.PropertyField(targetCamera, new GUIContent("Current Target Camera"));
				GUI.enabled = true;
			}
		}

		// Crash workaround
		for (int i = 0; i < ccmRenderers.Length; i++)
		{
			if (ccmRenderers[i].setRenderMode != ccmRenderers[i].renderMode)
			{
				ccmRenderers[i].DisableRendering();
				ccmRenderers[i].ReleaseEverything(false, false);
				ccmRenderers[i].resetButton = true;
			}
		}



		// Displacement Parameters
		if (MultiEditingNoControlMeshMode(ccmRenderers) == true)
		{
			EditorGUILayout.Space();
			EditorGUILayout.Space();
			EditorGUILayout.LabelField("Vertex Displacement", EditorStyles.boldLabel);
			EditorGUILayout.PropertyField(displacementEnabled, new GUIContent("Enable"));
			if (MultiEditingAllDisplacementEnabled(ccmRenderers) == true)
			{
				EditorGUILayout.PropertyField(displacement, new GUIContent("Height Map"));
				EditorGUILayout.PropertyField(displacementAmplitude, new GUIContent("Amplitude"));
				EditorGUILayout.PropertyField(displacementOffset, new GUIContent("Offset"));
				PropertyFloatSlider(displacementBase, 0.0f, 1.0f, new GUIContent("Base"));
			}
		}


		// Debug Parameters // CUSTOM DEMO
		EditorGUILayout.Space();
		EditorGUILayout.LabelField("Debug Parameters", EditorStyles.boldLabel);
		EditorGUILayout.PropertyField(enableDebugWireframe, new GUIContent("Enable Debug Wireframe"));
		EditorGUILayout.PropertyField(whiteWireframe, new GUIContent("White Wireframe"));
		EditorGUILayout.PropertyField(useVertexWelding, new GUIContent("Weld vertices at loading"));
		EditorGUILayout.PropertyField(loadFromUnityGPU, new GUIContent("Load from Unity Mesh on GPU"));
		EditorGUILayout.PropertyField(bakeCCMFileButton, new GUIContent("Bake CCM File"));
		EditorGUILayout.PropertyField(loadFromCCMFile, new GUIContent("Load from CCM File"));

		for (int i = 0; i < ccmRenderers.Length; i++)
		{
			ccmRenderers[i].usedCatmullClarkMaterial = ccmRenderers[i].catmullClarkRenderMaterial;

			ccmRenderers[i].usedCatmullClarkMaterial.SetFloat("u_WhiteWireframe", ccmRenderers[i].whiteWireframe == true ? 1.0f : 0.0f);

			if (ccmRenderers[i].enableDebugWireframe == true)
				ccmRenderers[i].usedCatmullClarkMaterial.EnableKeyword("CCM_WIREFRAME");
			else
				ccmRenderers[i].usedCatmullClarkMaterial.DisableKeyword("CCM_WIREFRAME");

			if (ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.ControlMesh == true)
				ccmRenderers[i].usedCatmullClarkMaterial.EnableKeyword("CCM_CONTROL_MODE");
			else
				ccmRenderers[i].usedCatmullClarkMaterial.DisableKeyword("CCM_CONTROL_MODE");

			if (ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.UniformTessellation == true)
				ccmRenderers[i].usedCatmullClarkMaterial.EnableKeyword("CCM_UNIFORM_MODE");
			else
				ccmRenderers[i].usedCatmullClarkMaterial.DisableKeyword("CCM_UNIFORM_MODE");
		}

		serializedObject.ApplyModifiedProperties();
	}



	private void PropertyIntSlider(SerializedProperty property, int leftValue, int rightValue, GUIContent label)
	{
		Rect rect = EditorGUILayout.GetControlRect();
		label = EditorGUI.BeginProperty(rect, label, property);

		EditorGUI.BeginChangeCheck();
		int newValue = EditorGUI.IntSlider(rect, label, property.intValue, leftValue, rightValue);
		// Only assign the value back if it was actually changed by the user.
		// Otherwise a single value will be assigned to all objects when multi-object editing,
		// even when the user didn't touch the control.
		if (EditorGUI.EndChangeCheck())
		{
			property.intValue = newValue;
		}
		EditorGUI.EndProperty();
		EditorGUILayout.Space();
	}

	private void PropertyFloatSlider(SerializedProperty property, float leftValue, float rightValue, GUIContent label)
	{
		Rect rect = EditorGUILayout.GetControlRect();
		label = EditorGUI.BeginProperty(rect, label, property);

		EditorGUI.BeginChangeCheck();
		float newValue = EditorGUI.Slider(rect, label, property.floatValue, leftValue, rightValue);
		// Only assign the value back if it was actually changed by the user.
		// Otherwise a single value will be assigned to all objects when multi-object editing,
		// even when the user didn't touch the control.
		if (EditorGUI.EndChangeCheck())
		{
			property.floatValue = newValue;
		}
		EditorGUI.EndProperty();
		EditorGUILayout.Space();
	}



	private int MultiEditingMinimumMaxDepth(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		int minimumMaxDepth = ccmRenderers[0].maxPossibleDepth;
		for (int i = 1; i < ccmRenderers.Length; i++)
			minimumMaxDepth = (int)Mathf.Min(ccmRenderers[i].maxPossibleDepth, minimumMaxDepth);
		return minimumMaxDepth;
	}

	private int MultiEditingMinimumSubdivisionDepth(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		int minimumMaxDepth = ccmRenderers[0].maxSubdivisionDepth;
		for (int i = 1; i < ccmRenderers.Length; i++)
			minimumMaxDepth = (int)Mathf.Min(ccmRenderers[i].maxSubdivisionDepth, minimumMaxDepth);
		return minimumMaxDepth;
	}

	private bool MultiEditingAllSameRenderMode(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		CatmullClarkMeshRenderer.RenderMode temp = ccmRenderers[0].renderMode;
		for (int i = 1; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].renderMode != temp)
				return false;
		return true;
	}

	private bool MultiEditingAllSameAutoLODMode(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		bool temp = ccmRenderers[0].uniformAutoLOD;
		for (int i = 1; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].uniformAutoLOD != temp)
				return false;
		return true;
	}

	private bool MultiEditingAllSetRenderMaterials(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		for (int i = 0; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].catmullClarkRenderMaterial == null)
				return false;
		return true;
	}

	private bool MultiEditingAllAdaptive(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		for (int i = 0; i < ccmRenderers.Length; i++)
			if (!(ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.AdaptiveTessellation || (ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.UniformTessellation && ccmRenderers[i].uniformAutoLOD == true)))
				return false;
		return true;
	}

	private bool MultiEditingNoControlMeshMode(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		for (int i = 0; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].renderMode == CatmullClarkMeshRenderer.RenderMode.ControlMesh)
				return false;
		return true;
	}

	private bool MultiEditingOneWireframeMode(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		for (int i = 0; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].debugMode == CatmullClarkMeshRenderer.ControlmeshDebugMode.Wireframe)
				return true;
		return false;
	}

	private string MultiEditingTotalMemoryConsumption(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		ulong totalByteSize = 0;
		for (int i = 0; i < ccmRenderers.Length; i++)
		{
			totalByteSize += ccmRenderers[i].byteSize;
		}

		if (false && totalByteSize >= (1 << 30))
		{
			return string.Format("{0} GB", totalByteSize >> 30);
		}
		else if (totalByteSize >= (1 << 20))
		{
			return string.Format("{0} MB", totalByteSize >> 20);
		}
		else if (totalByteSize >= (1 << 10))
		{
			return string.Format("{0} KB", totalByteSize >> 10);
		}
		else
		{
			return string.Format("{0} Bytes", totalByteSize);
		}
	}

	private int MultiEditingTotalTriangleCount(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		int totalTriangleCount = 0;
		for (int i = 0; i < ccmRenderers.Length; i++)
			totalTriangleCount += ccmRenderers[i].currentTriangleCount;
		return totalTriangleCount;
	}

	private bool MultiEditingAllDisplacementEnabled(CatmullClarkMeshRenderer[] ccmRenderers)
	{
		for (int i = 0; i < ccmRenderers.Length; i++)
			if (ccmRenderers[i].displacementEnabled == false)
				return false;
		return true;
	}
}
