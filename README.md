# neb3map-viewer
Drakensang Online Model and Map Viewer

![Preview](https://i.imgur.com/ETVwYrU.png)

A Nebula3-based model viewer using **GLFW3**, **GLAD**, and **DirectXTex** for legacy `.dds` decoding. \
Parses `.n3` and `.nvx2` assets — meshes from NVX2, materials and textures from N3.  \
As of **2025-10-16**, supports basic model viewing.  \
Full feature set (sfx + music, animations, particles, terrain [.map: decals, unlit, environment, water shaders], **FBX/Blender export**)? **Pay me.** \
Roadmap: postprocessing, pathfinding

Uses **deferred rendering**

### WORK DIR
> ⚠️ Project expects assets under:  
> `C:/drasa_online/work/`

### 📂 Supported / documented formats
**Animations:** `NAX3`, `NAC0`, `NA01`, `NA1C`  
**Meshes:** `NVX2`  
**Models:** `N3`  
**Maps:** `MAP`

### 🔧 Build Info
**Toolset:** v145  
**Language Standard:** C++17  
**Architecture:** x64  
**WIN SDK:** 10.0.26100.0
**Runtime:** MD (Release builds only)

### ❤️ Credits
**gyoerkaa** – NVX2 loader\
**Gscept** – Nebula3 source and references  