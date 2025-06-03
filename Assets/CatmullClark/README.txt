CatmullClarkMeshRenderer



DEPEDENCIES:
- Unity.Mathematics (Package Manager / Plus icon / Add from git URL / com.unity.mathematics)
- ShaderGraph 12.0 minimum



FEATURES:
- Unity implementation of "A Halfedge Refinement Rule for Parallel Catmull-Clark Subdivision" [J. Dupuy and K. Vanhoey, HPG 2021]
- Add-On component to Unity Renderer components that replaces their rendering with a Catmull Clark subdivision of the original mesh.
- Fast GPU computing of catmull clark subdivisions thanks to a halfedge internal representation.
- Adaptive tessellation through the pre-computed catmull clark subdivision levels using a Concurrent Binary Tree.
- GPU Skinning support for animated renderers
- Shadergraph CustomFunction implementation for compatibility accross render pipelines and versions



USAGE:
- Enable Read/Write on the source mesh files you want to use
- Create a new material using "Shader Graphs/CatmullClarkLit"
- Add the CatmullClarkMeshRenderer script to an object which has a Renderer already (either a MeshRenderer or a SkinnedMeshRenderer)
- Assign the new material to CatmullClarkMeshRenderer field, then enable the CatmullClarkMeshRenderer component



DETAILED USAGE:
- An enabled CatmullClarkMeshRenderer automatically replaces the UnityEngine.Renderer rendering through forceRenderingOff, do not disable the Unity Renderer component (it needs to be enabled for skinned animation updates)
- To create a new ShaderGraph for use with CatmullClarkMeshRenderers, copy the CatmullClarkVert CustomFunction node with its VertexID input node from the provided Assets/CatmullClark/Shaders/CatmullClarkLit.shadergraph and assign its output nodes to the Vertex stage of your new ShaderGraph. A Custom Interpolator needs to be added there to connect the "uv" output to. This Custom Interpolator needs to be used for any texture sampling UV input instead of the regular UV input method.
- Performance metrics either show the VRAM usage and current displayed triangle count for a single selected object, or the total counts for a selection of objects
- When either Adaptive Tessellation or Uniform Tessellation with auto LOD is used, the level of detail is automatically updated for the current camera (the scene view camera while editing, or the user set camera when in Play mode)
- Enabling the "Keep Quads" option of the Model Importer for your source mesh is encouraged to avoid changing the subdivision surface designed by the artist, as CatmullClarkMeshRenderer components support it.
- When the source mesh is modified through mesh editing or Model Importer settings, CatmullClarkMeshRenderer components using it need to be reset to reflect the change, either through disabling/enabling or changing the Rendering Mode parameter.



PARAMETERS:
- Force Editor Constant Updates - in Edit mode for the scene view camera, the LOD cannot update in real-time as MonoBehaviour updates happen sporadically, on middle or right mouse button up events in the scene view. Enabling this parameter solves that by forcing real-time constant editor MonoBehaviour.Update calls. It uses OnDrawGizmos() to schedule updates, so this only works when gizmos are enabled.

- Catmull Clark Material - the material to use for rendering with catmull clark subdivision. Must either use the "Shader Graphs/CatmullClarkLit" shader or a shadergraph that uses the required CatmullClarkVert CustomFunction.

- Subdivisions - how many catmull clark subdivision levels to compute. Maximum possible value is limited by the maximum possible buffer size allocation on the GPU.

- Rendering Mode - Adaptive Tessellation, Uniform Tessellation or Control Mesh

	- Adaptive Tessellation offers fully adaptive triangulation of the pre-computed catmull clark subdivision levels
		- Pause Adaptive Tessellation - freezes the adaptive triangulation update
		- Target Triangle Pixel Length - the triangle edge length in pixels to aim for while updating the adaptive triangulation

	- Uniform Tessellation renders the object with a specific depth level from the pre-computed catmull clark subdivision
		- Automatic LOD Selection - automatically switch between subdivision levels based on the object's bounds size in screen space
			-	LOD Bias - manually bias the automatic LOD selection
		- Level of Detail - when Automatic LOD Selection is disabled, sets which catmull clark subdivision level to display manually. There are always twice more levels available than the "Subvidisions" parameter value.

	- Control Mesh displays the control cage of the catmull clark subdivision with both triangles and quads and displays boundary edges.
		- Control Mesh Debug Mode - various debug display modes

- Play Mode Target Camera - when using either Adaptive Tessellation or Uniform Tessellation with Automatic LOD Selection, the level of detail is adaptive relative to this Camera while in Play mode (the object is still rendered for every camera).

- Current Target Camera (read-only) - displays the Camera currently used to update the adaptive level of detail. This is automatically set to the Scene view's camera while Editing, and switchted to the user specified Play Mode Target Camera when entering Play mode.

- Vertex Displacement - enable vertex displacement mapping using the HDRP Amplitude Parameterization
	- Height Map - the height map that you'd typically assign on the material directly for displacement mapping.
	- Amplitude - Set the amplitude of the Height Map.
	- Offset - Set the offset that HDRP applies to the Height Map.
	- Base - Use the slider to set the base for the Height Map.



LIMITATIONS:
- No semi-sharp creases as Unity's mesh loader doesn't support them. 
- Only one texture coordinate channel is supported. More could be added at a memory cost.
- Submeshes aren't supported: when using a CatmullClarkMeshRenderer, every submesh is merged into one halfedge mesh for rendering. Intended behaviour is merging the Quad and Triangle separate submeshes that are generated when the "Keep Quads" option is used, as both are supported simultaneously. Unintended consequence is having no support for separate submeshes with separate materials.



KNOWN ISSUES:
- Some meshes might be incompatible with subdivision. When this is the case, loading will fail with the "CatmullClarkMeshRenderer: Unsupported Geometry" error message.
- Currently there is no Per Object Motion vector generation => no TAA support.
- When used on an animated skinned mesh, there is currently a delay of a few frames in the displayed animation.
- CatmullClarkMeshRenderers are not rendered in Play mode when Paused.
- Objects rendered with CatmullClarkMeshRenderers are not selectable with a mouse click in the scene view. Disable the CatmullClarkMeshRenderer temporarily to do so.
- Rendering is bugged when using the experimental DirectX 12 Graphics API
