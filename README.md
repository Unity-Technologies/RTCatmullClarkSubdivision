# Catmull-Clark Mesh Renderer
This Unity Project is a demo of research work enabling parallellized Catmull-Clark subdivision on the GPU with adaptively-subdivised rendering.

It leverages the following works by [Jonathan Dupuy](https://onrendering.com/), [Thomas Deliot](https://thomasdeliot.wixsite.com/blog) and [Kenneth Vanhoey](http://www.kvanhoey.eu).
* _A Halfedge Refinement Rule for Parallel Catmull-Clark Subdivision_, published at High Performance Graphics 2021 ([paper](https://onrendering.com/data/papers/catmark/HalfedgeCatmullClark.pdf), [video](https://youtu.be/Ys2j5n9zWys), [presentation](https://www.youtube.com/watch?v=FFWgQZsfwy8))
* _Bisection Based Triangulation of Catmull Clark Subdivision_, Siggraph Talks 2022 ([paper](https://advances.realtimerendering.com/s2022/LebCatmullClark.pdf), [presentation](https://www.youtube.com/watch?v=p5tI2ddzMC4))

![Adaptive Mesh Subdivision!](/img/TRex-Unity.JPG "Adaptive TRex in Unity")

# Installation
In Unity Hub, add this folder and open with Unity 6 (tested with 6000.0.40f1 LTS and 6000.1.0a9_5daf0693406e_x64 ALPHA).

# Run Demo and Features
Open the `T-Rex Catmull-Clark` Asset.
By default, the mesh is subdivided on the GPU according to distance to camera: less refinement at a distance, more detailed for close-ups.

Options can be edited by
* Selecting `trex_walkcycle/Catmull_Clark_Trex_LOW` in the Hierarchy
* Go to the `Catmull Clark Mesh Renderer (Script)` tab in the Inspector

in which the following features are supported:
* Rendering Mode: Show Control Mesh (the coarsest level), Uniform or Adaptive Tessellation, respectively.
* Pause Adaptive Tessellation
* Set the desired Pixels per Triangle ratio

# Disclaimer
This demo software is provided "as is".
The repository owners cannot accept pull requests, GitHub review requests, or any other GitHub-hosted issue management requests.
