# neb3map-viewer
Drakensang Online Model & Map Viewer

In serious WIP!

The Nebula model viewew is simple glfw3 with glad loader,DirectXTex for old DDS parsing
The system is simple, it parses n3 and nvx2, takes the mesh from nvx2 and extracts the textures from n3, nothing else
all other FOURCC tags are skipped so the file can reads

Warning! The project work under C:/drasa_online/work/...

The nax2:
It parses the binary header, decodes the component mask, and extracts positions, normals, tangents, binormals, UVs, colors, and indices into a unified ObjVertex struct.
stores only geometry data! 

The n3:
linker for the nvx2 files, sets all the shaders, materials, animations etc
Handles enormous amount of FourCC data tags, skipping unneeded so it can gracefully load the model

FourCC structure of the .map binary:
The .map is custom Maya exported scene map with Radon Labs own plugin.
Thats are as follows: MOSD, IPAM, TRTS, TTES, TEVE, LAME, LPMT, TSNI, SNII, PORG, BVAN

Made with V143 toolset, C++17 for x64 
Thanks to: gyoerkaa for the nvx2loader! And Gscept!