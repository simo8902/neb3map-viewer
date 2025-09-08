/*                                         PRODBYSIMO
                          MADE WITH TOOLKIT V143 /MDd WIN 10.0.26100.0 SDK+

                                            09/01/2025
    * ================================ VERSION 0.0.2 START =============================
    * Finally made the parser respect every node in the .n3 instead of 
	* just grabbing the first mesh lol
    * Hooked up local -> world transform chains, so objects don't all chill at 0,0,0
    * Added primitive group support, so submeshes can render with the 
	* right material instead of being mashed togfether
    * Materials are now actually pulled from the node’s STXT tags, 
	* with fallbacks if they’re missing (white, black, flat normal)
    * Implemented UV transforms + UV set switching (yep, STUS/SCVU), 
	* pushed as uniforms into the shader
    * Shader updated with alpha cutout discard, so trees and 
	* leaves don’t look like cardboard anymore
    * Draw calls now go through a per-node DrawCmd, so each piece has its own mesh,
	* textures, transforms, and group index
    * Added a bunch of new tag handlers (joints, anim stuff, uv/scaling junk, etc.) 
	* so files stop spitting errors at random
    * Cleaned up & fixed some broken skips where the parser used to desync 
	* and break the whole model
    * ================================ VERSION 0.0.2 END ===============================
    
                                           09/02/2025
    * ================================ VERSION 0.0.3 START  ============================
    * added glm and nlohmann::json libs
    * TIRED ASF
    * ================================ VERSION 0.0.3 END ==============================='
    
    * ================================ VERSION 0.0.4 START  ============================
    * fixed dds loading
    * fixed some stuff
    * added more mess
    * json parser obviously (until I dont figure it out, I wont implement it)
    * fixed grid
    * TIRED ASF
    * ================================ VERSION 0.0.4 END ==============================='


    * ================================= TAGS TO PROCEED =================================
    * FX -> QRFE
    *
    * ================================= TAGS TO PROCEED END =============================
*/


/*
    BAKED TRANSFORM'S DATA
    =====================================================================================================
    ground_02 -> pre baked transform NVX2 and N3, dont resize! dont expect appearing as Transform in .map
    =====================================================================================================
*/

#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#include "json.hpp"
using json = nlohmann::json;

static const std::string MODELS_ROOT = "C:/drasa_online/work/models/";
static const std::string MESHES_ROOT = "C:/drasa_online/work/meshes/";

#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#undef APIENTRY
#include <DirectXTex/DirectXTex.h>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <direct.h>
#include <filesystem>
#include <array>
#include <random>
#include <unordered_set>
#include <sstream>

// LOG STRUCTURES
struct cerr_t {
    template<typename T>
    cerr_t& operator<<(const T& value) {
        std::cerr << "\033[38;2;255;0;0m" << value << "\033[0m";
        return *this;
    }

    cerr_t& operator<<(std::ostream& (*manip)(std::ostream&)) {
        std::cerr << manip;
        return *this;
    }
};

struct logout_t {
    template<typename shit>
    logout_t& operator<<(const shit& v) {
        std::cout << "\033[38;2;255;255;0m" << v << "\033[0m";
        return *this;
    }

    logout_t& operator<<(std::ostream& (*m)(std::ostream&)) {
        std::cerr << m;
        return *this;
    }
};

inline cerr_t errout;
inline logout_t logout;

std::unordered_map<std::string, GLuint> gTexCache;
GLuint gWhiteTex = 0;
GLuint gBlackTex = 0;
GLuint gFlatNormalTex = 0;

void ensureFallbacks();

GLuint loadNebulaByKeys(
    const std::unordered_map<std::string, std::string>& m,
    std::initializer_list<const char*> keys);

template<typename T>
T clamp(T v, T lo, T hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

GLuint gridVAO = 0, gridVBO = 0;
GLsizei gridVertexCount = 0;
constexpr double M_PI = 3.14159265358979323846;
unsigned int gDDSQuadVAO = 0, gDDSQuadVBO = 0, gDDSTexture = 0;
GLuint gDDSShader = 0;
bool gShowDDS = true;
bool gCPrev = false;

struct ObjVertex { float px, py, pz; float nx, ny, nz; float tx, ty, tz; float bx, by, bz; float u0, v0; float u1, v1; float cr, cg, cb, ca; float w0, w1, w2, w3; uint8_t j0, j1, j2, j3; };

struct Nvx2Group {
    uint32_t firstVertex = 0;
    uint32_t numVertices = 0;
    uint32_t firstTriangle = 0;
    uint32_t numTriangles = 0;
    uint32_t firstEdge = 0;
    uint32_t numEdges = 0;

    uint32_t firstIndex()   const { return firstTriangle * 3; }
    uint32_t indexCount()   const { return numTriangles * 3; }
    uint32_t baseVertex()   const { return firstVertex; }
};

struct Mesh {
    std::vector<ObjVertex> verts;
    std::vector<uint32_t> idx;
    std::vector<Nvx2Group> groups;
    GLuint vao = 0, vbo = 0, ebo = 0;
};

void setIdentity(float* mat) {
    mat[0] = 1.0f; mat[4] = 0.0f; mat[8] = 0.0f; mat[12] = 0.0f;
    mat[1] = 0.0f; mat[5] = 1.0f; mat[9] = 0.0f; mat[13] = 0.0f;
    mat[2] = 0.0f; mat[6] = 0.0f; mat[10] = 1.0f; mat[14] = 0.0f;
    mat[3] = 0.0f; mat[7] = 0.0f; mat[11] = 0.0f; mat[15] = 1.0f;
}

struct DrawCmd {
    Mesh mesh;
    std::string name;
    GLuint tex[4];
    bool   has[4];
    float  worldMatrix[16];
    int    group = -1;

    float uvXform[4][4];
    int   uvSet[4];

    void applyUVTransforms(GLuint shaderProgram) {
        for (int i = 0; i < 4; ++i) {
            std::string uniformName = "UvXform[" + std::to_string(i) + "]";
            GLint loc = glGetUniformLocation(shaderProgram, uniformName.c_str());
            if (loc >= 0) {
                glUniform4fv(loc, 1, uvXform[i]);
            }
        }
        GLint uvSetLoc = glGetUniformLocation(shaderProgram, "UvSet[0]");
        if (uvSetLoc >= 0) {
            glUniform1iv(uvSetLoc, 4, uvSet);
        }
    }

    DrawCmd() {
        memset(tex, 0, sizeof(tex));
        memset(has, 0, sizeof(has));
        setIdentity(worldMatrix);
        for (int s = 0; s < 4; ++s) {
            uvXform[s][0] = 1.0f;
            uvXform[s][1] = 1.0f;
            uvXform[s][2] = 0.0f;
            uvXform[s][3] = 0.0f;
            uvSet[s] = 0;
        }
    }
};

std::vector<DrawCmd> fromNodes;
std::vector<DrawCmd> withTransform;

float targetX = 0.0f, targetY = 0.0f, targetZ = 0.0f;

bool firstMouse = true;
double lastX = 400.0;
double lastY = 300.0;
bool mousePressed = false;
bool keys[1024];
float yaw = 0.0f;
float pitch = 0.0f;
float distance = 5.0f;
bool rightMouseDown = false;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

namespace fs = std::filesystem;

const std::vector<int> SUPPORTED_VERSIONS = { 1, 2 };

static void CheckShader(GLuint sh, const char* name) {
    GLint ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len = 0; glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0'); glGetShaderInfoLog(sh, len, nullptr, log.data());
        errout << "[SHADER COMPILE ERROR] " << name << "\n" << log << "\n";
    }
}
static void CheckProgram(GLuint prog) {
    GLint ok = 0; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0; glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0'); glGetProgramInfoLog(prog, len, nullptr, log.data());
        errout << "[PROGRAM LINK ERROR]\n" << log << "\n";
    }
}


std::string resolveTexture(const std::string& texName) {
    std::string clean = texName;
    if (clean.rfind("tex:", 0) == 0) clean = clean.substr(4);
    return "C:/drasa_online/work/textures/" + clean + ".dds";
}
#ifndef GL_COMPRESSED_RGBA_S3TC_DXT1_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT1_EXT 0x83F1
#endif
#ifndef GL_COMPRESSED_RGBA_S3TC_DXT3_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT3_EXT 0x83F2
#endif
#ifndef GL_COMPRESSED_RGBA_S3TC_DXT5_EXT
#define GL_COMPRESSED_RGBA_S3TC_DXT5_EXT 0x83F3
#endif
#ifndef GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT
#define GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT 0x8C4D
#endif
#ifndef GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT
#define GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT 0x8C4E
#endif
#ifndef GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT
#define GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT 0x8C4F
#endif

bool debug_output = true;

GLuint LoadDDS(const std::string& path) {
    auto itc = gTexCache.find(path);
    if (itc != gTexCache.end()) return itc->second;

    std::string norm = path;
    for (auto& ch : norm) if (ch == '/') ch = '\\';

    {
        std::ifstream probe(norm, std::ios::binary);
        if (!probe.good()) {
            if (debug_output) errout << "[TEX ERR] Not found / unreadable: " << norm << "\n";
            return 0;
        }
    }

    DirectX::ScratchImage image;
    DirectX::TexMetadata meta;
    HRESULT hr = E_FAIL;

    hr = DirectX::LoadFromDDSFile(
        std::wstring(norm.begin(), norm.end()).c_str(),
        DirectX::DDS_FLAGS_NONE,
        &meta, image
    );

    if (FAILED(hr)) {
        hr = DirectX::LoadFromDDSFile(
            std::wstring(norm.begin(), norm.end()).c_str(),
            DirectX::DDS_FLAGS_LEGACY_DWORD | DirectX::DDS_FLAGS_NO_LEGACY_EXPANSION,
            &meta, image
        );
    }

    if (FAILED(hr)) {
        std::ifstream f(norm, std::ios::binary | std::ios::ate);
        if (!f) {
            if (debug_output) errout << "[TEX ERR] Failed to open for memory load: " << norm << "\n";
            return 0;
        }
        std::streamsize sz = f.tellg();
        f.seekg(0, std::ios::beg);
        std::vector<uint8_t> bytes(sz);
        if (!f.read(reinterpret_cast<char*>(bytes.data()), sz)) {
            if (debug_output) errout << "[TEX ERR] Failed to read bytes for memory load: " << norm << "\n";
            return 0;
        }
        hr = DirectX::LoadFromDDSMemory(bytes.data(), bytes.size(),
            DirectX::DDS_FLAGS_NONE, &meta, image);
        if (FAILED(hr)) {
            hr = DirectX::LoadFromDDSMemory(bytes.data(), bytes.size(),
                DirectX::DDS_FLAGS_LEGACY_DWORD | DirectX::DDS_FLAGS_NO_LEGACY_EXPANSION,
                &meta, image);
        }
    }

    if (FAILED(hr)) {
        if (debug_output) errout << "[TEX ERR] Failed to load DDS: " << norm << "\n";
        return 0;
    }

    const DirectX::Image* img = image.GetImage(0, 0, 0);
    if (!img) {
        if (debug_output) errout << "[TEX ERR] No image data: " << norm << "\n";
        return 0;
    }

    GLenum glFormat = 0;
    bool isSRGB = false;

    switch (meta.format) {
    case DXGI_FORMAT_BC1_UNORM:        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT; break;
    case DXGI_FORMAT_BC1_UNORM_SRGB:   glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT; isSRGB = true; break;
    case DXGI_FORMAT_BC2_UNORM:        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT3_EXT; break;
    case DXGI_FORMAT_BC2_UNORM_SRGB:   glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT; isSRGB = true; break;
    case DXGI_FORMAT_BC3_UNORM:        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT5_EXT; break;
    case DXGI_FORMAT_BC3_UNORM_SRGB:   glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT; isSRGB = true; break;

    case DXGI_FORMAT_R8G8B8A8_UNORM:
    case DXGI_FORMAT_R8G8B8A8_UNORM_SRGB:
    default:
    {
        DirectX::ScratchImage conv;
        const DirectX::Image* src = img;
        DXGI_FORMAT want = DXGI_FORMAT_R8G8B8A8_UNORM;
        if (meta.format == DXGI_FORMAT_R8G8B8A8_UNORM_SRGB) isSRGB = true;

        if (meta.format != DXGI_FORMAT_R8G8B8A8_UNORM &&
            meta.format != DXGI_FORMAT_R8G8B8A8_UNORM_SRGB) {
            HRESULT hr2 = DirectX::Decompress(*img, DXGI_FORMAT_R8G8B8A8_UNORM, conv);
            if (FAILED(hr2)) {
                if (debug_output) errout << "[TEX ERR] Decompress failed: " << norm << "\n";
                return 0;
            }
            src = conv.GetImage(0, 0, 0);
        }

        GLuint tex = 0;
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_2D, tex);
        glTexImage2D(GL_TEXTURE_2D, 0, isSRGB ? GL_SRGB8_ALPHA8 : GL_RGBA8,
            (GLsizei)src->width, (GLsizei)src->height,
            0, GL_RGBA, GL_UNSIGNED_BYTE, src->pixels);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glBindTexture(GL_TEXTURE_2D, 0);

        gTexCache[path] = tex;
        if (debug_output) errout << "[TEX] Loaded RGBA8: " << norm << " -> " << tex << "\n";
        return tex;
    }
    }

    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    for (size_t level = 0; level < meta.mipLevels; ++level) {
        const DirectX::Image* mip = image.GetImage(level, 0, 0);
        glCompressedTexImage2D(GL_TEXTURE_2D, (GLint)level, glFormat,
            (GLsizei)mip->width, (GLsizei)mip->height, 0,
            (GLsizei)mip->slicePitch, mip->pixels);
    }
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        meta.mipLevels > 1 ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    bool isFaceLike = (norm.find("mask") != std::string::npos || norm.find("face") != std::string::npos);
    GLenum wrapMode = isFaceLike ? GL_CLAMP_TO_EDGE : GL_REPEAT;
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrapMode);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrapMode);

    glBindTexture(GL_TEXTURE_2D, 0);

    gTexCache[path] = tex;
    if (debug_output) std::cout << "[TEX] Loaded BCn: " << norm << " -> " << tex << "\n";
    return tex;
}



struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
};

struct Vec4 {
    float x, y, z, w;
    Vec4() : x(0), y(0), z(0), w(1) {}
    Vec4(float x_, float y_, float z_, float w_) : x(x_), y(y_), z(z_), w(w_) {}
};


void multiplyMatrices(float* r, const float* a, const float* b) {
    for (int row = 0; row < 4; ++row) {
        for (int col = 0; col < 4; ++col) {
            r[col * 4 + row] =
                a[0 * 4 + row] * b[col * 4 + 0] +
                a[1 * 4 + row] * b[col * 4 + 1] +
                a[2 * 4 + row] * b[col * 4 + 2] +
                a[3 * 4 + row] * b[col * 4 + 3];
        }
    }
}

void translateMatrix(float* matrix, float x, float y, float z) {
    float translation[16] = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f,
        x, y, z, 1.0f
    };
    float temp[16];
    multiplyMatrices(temp, matrix, translation);
    for (int i = 0; i < 16; i++) matrix[i] = temp[i];
}

void rotateYMatrix(float* matrix, float angle) {
    float cosA = cos(angle);
    float sinA = sin(angle);
    float rotation[16] = {
        cosA, 0.0f, -sinA, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f,
        sinA, 0.0f, cosA, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    float temp[16];
    multiplyMatrices(temp, matrix, rotation);
    for (int i = 0; i < 16; i++) matrix[i] = temp[i];
}

void rotateXMatrix(float* matrix, float angle) {
    float cosA = cos(angle);
    float sinA = sin(angle);
    float rotation[16] = {
        1.0f, 0.0f, 0.0f, 0.0f,
        0.0f, cosA, sinA, 0.0f,
        0.0f, -sinA, cosA, 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f
    };
    float temp[16];
    multiplyMatrices(temp, matrix, rotation);
    for (int i = 0; i < 16; i++) matrix[i] = temp[i];
}

void createPerspectiveMatrix(float* m, float fov, float aspect, float znear, float zfar) {
    float f = 1.0f / tanf(fov * static_cast<float>(M_PI) / 180.0f / 2.0f);
    for (int i = 0; i < 16; i++) m[i] = 0;
    m[0] = f / aspect;
    m[5] = f;
    m[10] = (zfar + znear) / (znear - zfar);
    m[11] = -1;
    m[14] = (2 * zfar * znear) / (znear - zfar);
}


void createLookAtMatrix(float* m,
    float eyeX, float eyeY, float eyeZ,
    float cx, float cy, float cz,
    float upX, float upY, float upZ) {

    float f[3] = { cx - eyeX, cy - eyeY, cz - eyeZ };
    float fl = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    f[0] /= fl; f[1] /= fl; f[2] /= fl;

    float up[3] = { upX, upY, upZ };
    float s[3] = {
        f[1] * up[2] - f[2] * up[1],
        f[2] * up[0] - f[0] * up[2],
        f[0] * up[1] - f[1] * up[0]
    };
    float sl = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    s[0] /= sl; s[1] /= sl; s[2] /= sl;

    float u[3] = {
        s[1] * f[2] - s[2] * f[1],
        s[2] * f[0] - s[0] * f[2],
        s[0] * f[1] - s[1] * f[0]
    };

    m[0] = s[0]; m[4] = s[1]; m[8] = s[2];  m[12] = -(s[0] * eyeX + s[1] * eyeY + s[2] * eyeZ);
    m[1] = u[0]; m[5] = u[1]; m[9] = u[2];  m[13] = -(u[0] * eyeX + u[1] * eyeY + u[2] * eyeZ);
    m[2] = -f[0]; m[6] = -f[1]; m[10] = -f[2]; m[14] = f[0] * eyeX + f[1] * eyeY + f[2] * eyeZ;
    m[3] = 0;    m[7] = 0;    m[11] = 0;    m[15] = 1;
}

void quatToMatrix(float* mat, const Vec4& q) {
    float x = q.x, y = q.y, z = q.z, w = q.w;
    float x2 = x + x, y2 = y + y, z2 = z + z;
    float xx = x * x2, xy = x * y2, xz = x * z2;
    float yy = y * y2, yz = y * z2, zz = z * z2;
    float wx = w * x2, wy = w * y2, wz = w * z2;

    mat[0] = 1.0f - (yy + zz); mat[4] = xy - wz;         mat[8] = xz + wy;         mat[12] = 0.0f;
    mat[1] = xy + wz;          mat[5] = 1.0f - (xx + zz); mat[9] = yz - wx;         mat[13] = 0.0f;
    mat[2] = xz - wy;          mat[6] = yz + wx;         mat[10] = 1.0f - (xx + yy); mat[14] = 0.0f;
    mat[3] = 0.0f;             mat[7] = 0.0f;            mat[11] = 0.0f;            mat[15] = 1.0f;
}

void buildTransformMatrix(float* m, const Vec3& p, const Vec4& q, const Vec3& s) {
    float S[16], R[16], T[16], tmp[16];
    setIdentity(S); setIdentity(R); setIdentity(T);
    S[0] = s.x; S[5] = s.y; S[10] = s.z;
    quatToMatrix(R, q);
    T[12] = p.x; T[13] = p.y; T[14] = p.z;
    multiplyMatrices(tmp, R, S);
    multiplyMatrices(m, T, tmp);
}

struct Node {
    std::string node_name;
    std::string node_type;
    std::shared_ptr<Node> node_parent;
    std::vector<std::shared_ptr<Node>> node_children;

    std::unordered_map<std::string, std::string> shader_textures;
    std::unordered_map<std::string, float> shader_parameters;
    std::unordered_map<int, std::array<float, 4>> uvXformBySlot;
    std::unordered_map<int, int> uvSetBySlot;

    std::string mesh_resource_id;
    std::string model_node_type;
    int primitive_group_idx;
    std::unordered_map<int, std::vector<int>> skin_fragments;

    Vec3 position;
    Vec3 scale;
    Vec4 rotation;

    float localMatrix[16];
    float worldMatrix[16];
    bool worldMatrixDirty;

    Node() : primitive_group_idx(-1), position(0, 0, 0), scale(1, 1, 1), rotation(0, 0, 0, 1), worldMatrixDirty(true) {
        setIdentity(localMatrix);
        setIdentity(worldMatrix);
    }
    void buildLocalMatrix() {
        buildTransformMatrix(localMatrix, position, rotation, scale);
        worldMatrixDirty = true;
    }

    void computeWorldMatrix() {
        if (!worldMatrixDirty) return;
        if (node_parent) {
            node_parent->computeWorldMatrix();
            float tmp[16];
            multiplyMatrices(tmp, node_parent->worldMatrix, localMatrix);
            for (int i = 0; i < 16; ++i) worldMatrix[i] = tmp[i];
        }
        else {
            for (int i = 0; i < 16; ++i) worldMatrix[i] = localMatrix[i];
        }
        worldMatrixDirty = false;
    }
};

void buildTransformMatrixColumnMajor(float* mat, const Vec3& pos, const Vec4& rot, const Vec3& scale) {
    float x = rot.x, y = rot.y, z = rot.z, w = rot.w;
    float x2 = x + x, y2 = y + y, z2 = z + z;
    float xx = x * x2, xy = x * y2, xz = x * z2;
    float yy = y * y2, yz = y * z2, zz = z * z2;
    float wx = w * x2, wy = w * y2, wz = w * z2;

    mat[0] = (1.0f - (yy + zz)) * scale.x;
    mat[1] = (xy + wz) * scale.x;
    mat[2] = (xz - wy) * scale.x;
    mat[3] = 0.0f;

    mat[4] = (xy - wz) * scale.y;
    mat[5] = (1.0f - (xx + zz)) * scale.y;
    mat[6] = (yz + wx) * scale.y;
    mat[7] = 0.0f;

    mat[8] = (xz + wy) * scale.z;
    mat[9] = (yz - wx) * scale.z;
    mat[10] = (1.0f - (xx + yy)) * scale.z;
    mat[11] = 0.0f;

    mat[12] = pos.x;
    mat[13] = pos.y;
    mat[14] = pos.z;
    mat[15] = 1.0f;

}



class Parser {
public:

    std::vector<std::shared_ptr<Node>> cleanNodes() const {
        std::vector<std::shared_ptr<Node>> out;
        std::unordered_set<std::string> seen;

        auto hashNode = [](const Node& n) {
            std::ostringstream oss;
            oss << n.mesh_resource_id << "#" << n.primitive_group_idx;
            for (int i = 0; i < 16; i++) oss << "," << n.worldMatrix[i];
            return oss.str();
            };

        for (auto& n : n3node_list) {
            if (n->mesh_resource_id.empty()) continue;
            if (!n->node_children.empty()) continue;

            std::string key = hashNode(*n);
            if (seen.insert(key).second) {
                out.push_back(n);
            }
        }
        return out;
    }

    template<typename T>
    T read_n3_value() {
        T v{};
        this->n3file.read(reinterpret_cast<char*>(&v), sizeof(T));
        return v;
    }
    static const std::unordered_set<std::string>& ValidFourCC() {
        static const std::unordered_set<std::string> s{
            "LDOM","CHRN","LKSN","LNKS","ANIM","TNJN","NJNU",
            "JONT","NOJA","NOJB","NOJC","NOJD","NOJE","NOJF",
            "NOJG","NOJH","NOJI","NOJJ","NOJK","NOJL","NOJM",
            "NOJN","NOJO","NOJP","NOJQ","NOJR","NOJS","NOJT",
            "TRAV","NFRT","LBOX","NSHC","TXTS","TLFS","CEVS",
            "HSAC","RDHS","PTNM","HSEM","SEMS","IRGP","FKSN",
            "GRFS", "TNIS","DNM<","DNM>","LDM>","LDM<","EOF_","_FOE"
        };
        return s;
    }

    bool should_reverse_fourcc(const std::string& raw) {
        static const std::unordered_set<std::string> reversed_tags{
            "3BEN", "MODL", "CHRN", "MINA", "TNOJ", "XOBL", "VART",
            "NTJN", "UNJN", "TRFN", "CHSN", "MESH", "PRIG", "NSKF",
            "SFRG", "SHDR", "MNTP", "SMES", "CASH"
        };
        return reversed_tags.count(raw) > 0;
    }


private:

    bool is_valid_fourcc(const std::string& t) {
        return t.size() == 4 && ValidFourCC().count(t) > 0;
    }



    std::string filepath;
    std::ifstream n3file;
    std::string byteorder;
    int n3version = 0;
    std::string n3modeltype;
    std::string n3modelname;
    std::vector<std::shared_ptr<Node>> n3node_list;

    std::string read_n3_string() {
        uint16_t len = read_n3_value<uint16_t>();
        std::string s;
        if (len) {
            s.resize(len);
            n3file.read(&s[0], len);
        }
        return s;
    }

    std::string read_n3_fourcc() {
        char b[4];
        if (!n3file.read(b, 4)) throw std::runtime_error("eof");
        std::string tag(b, 4);
        if (should_reverse_fourcc(tag)) {
            std::reverse(tag.begin(), tag.end());
        }
        return tag;
    }

    bool looks_like_tag(const std::string& t) {
        if (t.size() != 4) return false;
        for (char c : t) {
            if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '<' || c == '>'))
                return false;
        }
        return true;
    }

    bool peek_known_four(std::streampos at, std::string& out) {
        if (at < 0) return false;
        std::streampos p = n3file.tellg();
        n3file.clear();
        n3file.seekg(at);

        char b[4];
        if (!n3file.read(b, 4)) { n3file.clear(); n3file.seekg(p); return false; }

        std::string t(b, 4);
        n3file.clear();
        n3file.seekg(p);
        if (!is_valid_fourcc(t)) return false;
        out = std::move(t);
        return true;
    }
    bool resync_to_next_known(std::size_t maxScan = 1 << 20) {
        std::streampos start = n3file.tellg();
        for (std::size_t i = 0; i < maxScan && n3file.good(); ++i) {
            std::string t;
            std::streampos pos = start + static_cast<std::streamoff>(i);
            if (!peek_known_four(pos, t)) continue;
            n3file.clear();
            n3file.seekg(pos);
            return true;
        }
        n3file.clear();
        n3file.seekg(start);
        return false;
    }

    std::string peek_fourcc(std::ifstream& f, const std::string& byteorder) {
        std::streampos p = f.tellg();
        char b[4];
        if (!f.read(b, 4)) {
            f.clear();
            f.seekg(p);
            return "";
        }
        f.clear();
        f.seekg(p);

        return std::string(b, 4);
    }

    bool peek_four(std::streampos at, std::string& out) {
        if (at < 0) return false;
        std::streampos p = n3file.tellg();
        n3file.clear();
        n3file.seekg(at);

        char b[4];
        if (!n3file.read(b, 4)) {
            n3file.clear();
            n3file.seekg(p);
            return false;
        }
        out.assign(b, 4);
        n3file.clear();
        n3file.seekg(p);
        return true;
    }
    bool resync_to_valid(std::size_t maxScan = 1 << 20) {
        std::streampos start = n3file.tellg();
        for (std::size_t i = 0; i < maxScan && n3file.good(); ++i) {
            std::string t; std::streampos pos = start + static_cast<std::streamoff>(i);
            if (!peek_four(pos, t)) break;
            if (ValidFourCC().count(t)) {
                n3file.clear();
                n3file.seekg(pos);
                return true;
            }
        }
        n3file.clear();
        n3file.seekg(start);
        return false;
    }
    bool is_valid_here() {
        std::string t; std::streampos p = n3file.tellg();
        if (!peek_four(p, t)) return false;
        return ValidFourCC().count(t) != 0;
    }

    bool looks_allowed(const std::string& t) {
        if (t.size() != 4) return false;
        for (char c : t) if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '<' || c == '>' || c == '_')) return false;
        return true;
    }


    bool is_known_here() {
        std::string t; std::streampos p = n3file.tellg();
        if (!peek_four(p, t)) return false;
        return looks_allowed(t);
    }


    bool resync(std::size_t maxScan = 1 << 20) {
        std::streampos start = n3file.tellg();
        for (std::size_t i = 0; i < maxScan && n3file.good(); ++i) {
            std::string t; std::streampos pos = start + static_cast<std::streamoff>(i);
            if (!peek_four(pos, t)) break;
            if (looks_allowed(t)) { n3file.clear(); n3file.seekg(pos); return true; }
        }
        n3file.clear(); n3file.seekg(start); return false;
    }

    bool try_skip_bytes(std::size_t n) {
        std::streampos p = n3file.tellg();
        n3file.seekg(static_cast<std::streamoff>(n), std::ios::cur);
        if (!n3file.good()) { n3file.clear(); n3file.seekg(p); return false; }
        if (is_valid_here()) return true;
        n3file.clear(); n3file.seekg(p);
        return false;
    }


    bool skip_envelope() { for (int i = 0; i < 8; i++) (void)read_n3_value<float>(); (void)read_n3_value<int32_t>(); return true; }
    void skip_floats(int n) { for (int i = 0; i < n; i++) (void)read_n3_value<float>(); }
    void skip_ints(int n) { for (int i = 0; i < n; i++) (void)read_n3_value<int32_t>(); }

    bool generic_skip() {
        std::streampos p = n3file.tellg();
        if (p == std::streampos(-1)) return false;

        {
            std::streampos q = n3file.tellg();
            int32_t sz = read_n3_value<int32_t>();
            if (sz > 0 && sz < 10'000'000) {
                n3file.seekg(static_cast<std::streamoff>(sz), std::ios::cur);
                if (n3file.good() && is_valid_here()) return true;
            }
            n3file.clear(); n3file.seekg(q);
        }

        {
            std::streampos q = n3file.tellg();
            uint16_t L = read_n3_value<uint16_t>();
            if (L < 65535) {
                n3file.seekg(L, std::ios::cur);
                if (n3file.good() && is_valid_here()) return true;
            }
            n3file.clear(); n3file.seekg(q);
        }

        {
            std::streampos q = n3file.tellg();
            auto rd = [&](uint16_t& L) {
                L = read_n3_value<uint16_t>();
                n3file.seekg(L, std::ios::cur);
                };
            uint16_t L1 = 0, L2 = 0;
            rd(L1); rd(L2);
            if (n3file.good() && is_valid_here()) return true;
            n3file.clear(); n3file.seekg(q);
        }

        n3file.clear(); n3file.seekg(p);
        return resync_to_valid();
    }

public:

    int currentSlot = 0;

    Parser() {}
    const std::vector<std::shared_ptr<Node>>& getNodes() const { return n3node_list; }
    void report(const std::string& t, const std::string& m) { errout << "[" << t << "] " << m << std::endl; }

    int infer_slot_from_key(const std::string& k) {
        if (!k.empty() && isdigit((unsigned char)k.back())) return k.back() - '0';
        std::string low = k;
        std::transform(low.begin(), low.end(), low.begin(), [](unsigned char c) { return (char)tolower(c); });
        if (low.find("diff") != std::string::npos || low.find("albedo") != std::string::npos || low.find("color") != std::string::npos) return 0;
        if (low.find("spec") != std::string::npos) return 1;
        if (low.find("bump") != std::string::npos || low.find("normal") != std::string::npos) return 2;
        if (low.find("emis") != std::string::npos) return 3;
        return 0;
    }

    bool parse_tag_shape(const std::string& t, std::shared_ptr<Node>& node) {
        if (t == "MESH") { node->mesh_resource_id = read_n3_string(); std::cout << "        mesh_res_id=" << node->mesh_resource_id << std::endl; return true; }
        if (t == "PGRI") {
            node->primitive_group_idx = read_n3_value<int32_t>();
            std::cout << "        primitive_group_idx=" << node->primitive_group_idx << std::endl;
            return true;
        }
        return false;
    }

    void parse_NJNT_block() {
        int32_t num = read_n3_value<int32_t>();
        if (num < 0 || num > 100000) throw std::runtime_error("NJNT count");
        for (int32_t k = 0; k < num; ++k) {
            std::string nm = read_n3_string();
            std::string tag = read_n3_fourcc();
            if (tag != "TNOJ") throw std::runtime_error("TNOJ");
            uint32_t ji = read_n3_value<uint32_t>();
            int32_t  pi = read_n3_value<int32_t>();
            float t0 = read_n3_value<float>(); float t1 = read_n3_value<float>(); float t2 = read_n3_value<float>();
            float q0 = read_n3_value<float>(); float q1 = read_n3_value<float>(); float q2 = read_n3_value<float>(); float q3 = read_n3_value<float>();
            float s0 = read_n3_value<float>(); float s1 = read_n3_value<float>(); float s2 = read_n3_value<float>();
            (void)nm; (void)ji; (void)pi; (void)t0; (void)t1; (void)t2; (void)q0; (void)q1; (void)q2; (void)q3; (void)s0; (void)s1; (void)s2;
        }
    }


    // NEBULA DEVICE BACKEND
    // -- characterskin;characterskinnode
    bool skip4cc(const std::string& t, std::shared_ptr<Node>& node) {
        auto skipJoint = [&](int f) {
            (void)read_n3_value<int32_t>();
            (void)read_n3_value<int32_t>();
            for (int i = 0; i < f; i++) (void)read_n3_value<float>();
            (void)read_n3_string();
            return true;
            };

        if (t == "NSKL") {
            int32_t numSkinLists = read_n3_value<int32_t>();
            std::cout << "NSKL: numSkinLists=" << numSkinLists << std::endl;
            return true;
        }
        else if (t == "LNKS") {
            std::string skinListName = read_n3_string();
            int32_t numSkins = read_n3_value<int32_t>();
            std::cout << "LNKS: " << skinListName << " numSkins=" << numSkins << std::endl;
            for (int j = 0; j < numSkins; j++) {
                std::string skinName = read_n3_string();
                std::cout << "  skin[" << j << "]=" << skinName << std::endl;
            }
            std::streampos before = n3file.tellg();
            std::string next = peek_fourcc(n3file, byteorder);
            if (!looks_like_tag(next)) (void)read_n3_value<uint8_t>();
            else { n3file.clear(); n3file.seekg(before); }
            return true;
        }
        else if (t == "CHSN") {
            node->node_type = "CHSN";
            return true;
        }

      
        else if (t == "SKNL") {
            std::string skinListName = read_n3_string();
            int32_t numSkins = read_n3_value<int32_t>();
            std::cout << "SKNL: name=" << skinListName << " numSkins=" << numSkins << std::endl;

            for (int i = 0; i < numSkins && i < 1000; i++) { 
                std::string skinName = read_n3_string();
                std::cout << "  skin[" << i << "]: " << skinName << std::endl;
            }

            if (n3version >= 2) {
                std::streampos before = n3file.tellg();
                auto next = peek_fourcc(n3file, byteorder);
                if (!looks_like_tag(next)) {
                    (void)read_n3_value<uint8_t>();
                    (void)read_n3_string();
                }
                else {
                    n3file.clear();
                    n3file.seekg(before);
                }
            }
            return true;
        }
        else if (t == "NJNT") {
            parse_NJNT_block();
            return true;
        }
        else if (t == "JONT") {
            int32_t jointIndex = read_n3_value<int32_t>();
            int32_t parentJointIndex = read_n3_value<int32_t>();

            float poseTranslation[4];
            for (int i = 0; i < 4; i++) {
                poseTranslation[i] = read_n3_value<float>();
            }

            float poseRotation[4];
            for (int i = 0; i < 4; i++) {
                poseRotation[i] = read_n3_value<float>();
            }

            float poseScale[4];
            for (int i = 0; i < 4; i++) {
                poseScale[i] = read_n3_value<float>();
            }

            std::string jointName = read_n3_string();

            std::cout << "JONT: idx=" << jointIndex << " parent=" << parentJointIndex
                << " name=" << jointName << std::endl;
            return true;
        }

        else if (t == "ANIM" || t == "MINA" || t == "VART" || t == "TRAV") {
            std::string resourceId = read_n3_string();
            std::cout << t << ": " << resourceId << std::endl;
            return true;
        }
        else if (t == "NJMS" || t == "TNJN") {
            int32_t value = read_n3_value<int32_t>();
            std::cout << t << ": " << value << std::endl;
            return true;
        }
        else if (t == "BJNT") return skipJoint(0);
        else if (t == "SJNT") return skipJoint(11);
        else if (t == "TNOJ") {
            uint32_t idx = read_n3_value<uint32_t>();
            uint32_t parentIdx = read_n3_value<uint32_t>();
            float m[16];
            n3file.read(reinterpret_cast<char*>(m), sizeof(m));
            std::cout << "TNOJ bone idx=" << idx << " parent=" << parentIdx << std::endl;
            return true;
        }

        else if (t == "CASH") { (void)read_n3_value<uint8_t>(); return true; }
        else if (t == "SHDR") { (void)read_n3_string(); return true; }

        else if (t == "LBOX") { for (int i = 0; i < 8; i++) (void)read_n3_value<float>(); return true; }
        else if (t == "RPIV" || t == "VIPR" || t == "SPIV" || t == "VIPS") { for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }

        else if (t == "SVSP" || t == "PSVS" || t == "SLKV" || t == "VKLS" ||
            t == "SLOP" || t == "SROF" || t == "SBBO" || t == "SSTS" || t == "SVAF") {
            (void)read_n3_value<uint8_t>();
            return true;
        }

        else if (t == "SSTA" || t == "SLPT" || t == "ANNO" || t == "SPNM" || t == "SVCN" ||
            t == "SANI" || t == "MNTP" || t == "SMSH" || t == "SSHD" || t == "SVRT" || t == "SANM") {
            (void)read_n3_string();
            return true;
        }

        else if (t == "SINT") { (void)read_n3_string(); (void)read_n3_value<int32_t>(); return true; }
        else if (t == "SFLT") { (void)read_n3_string(); (void)read_n3_value<float>(); return true; }
        else if (t == "SVEC") { (void)read_n3_string(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }

        else if (t == "SFRG") {
            int group_idx = read_n3_value<int32_t>();
            int num = read_n3_value<int32_t>();
            std::vector<int> joints(num);
            for (int i = 0; i < num; ++i) joints[i] = read_n3_value<int32_t>();
            node->skin_fragments[group_idx] = std::move(joints);
            return true;
        }

        else if (t == "BASE" || t == "ESAB" || t == "SAGR" || t == "RGAS" ||
            t == "PGRI" || t == "NSKF" || t == "BGFR" || t == "SFGI" ||
            t == "BGJP" || t == "SJID") {
            (void)read_n3_value<int32_t>();
            return true;
        }

        else if (t == "SMID" || t == "DIMS" || t == "SMAD" || t == "DAMS" || t == "SEMD" ||
            t == "SACD" || t == "SRMN" || t == "SRMX" || t == "SGRV" || t == "SPST" ||
            t == "SCVR" || t == "SCVS" || t == "SCVT" || t == "SPCT" || t == "STDL") {
            (void)read_n3_value<float>();
            return true;
        }

        else if (t == "STUS") {
            int32_t slot = read_n3_value<int32_t>();
            float us = read_n3_value<float>();
            float vs = read_n3_value<float>();
            float uo = read_n3_value<float>();
            float vo = read_n3_value<float>();
            currentSlot = slot;
            node->uvXformBySlot[slot] = std::array<float, 4>{ us,vs,uo,vo };
            if (!node->uvSetBySlot.count(slot)) node->uvSetBySlot[slot] = 0;
            return true;
        }

        else if (t == "SCVU") {
            int32_t slot = read_n3_value<int32_t>();
            uint8_t uvset = read_n3_value<uint8_t>();
            (void)read_n3_value<uint8_t>();
            (void)read_n3_value<uint8_t>();
            (void)read_n3_value<uint8_t>();
            node->uvSetBySlot[slot] = (int)uvset;
            currentSlot = slot;
            return true;
        }

        else if (t == "SSPI") { (void)read_n3_value<int32_t>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }

        else if (t == "EFRQ" || t == "PLFT" || t == "PSMN" || t == "PSMX" || t == "PSVL" ||
            t == "PRVL" || t == "PSZE" || t == "PMSS" || t == "PTMN" || t == "PVLF" ||
            t == "PAIR" || t == "PRED" || t == "PGRN" || t == "PBLU" || t == "PALP" ||
            t.rfind("SCV", 0) == 0) {
            for (int i = 0; i < 8; i++) (void)read_n3_value<float>();
            (void)read_n3_value<int32_t>();
            return true;
        }

        else if (t == "LDOM") {
            for (int i = 0; i < 8; i++) (void)read_n3_value<float>();
            (void)read_n3_value<int32_t>();
            return true;
        }
        else if (t == "NRHC"  || t == "TNJN" || t == "NJNU" ||
            t == "NOJA" || t == "NOJB" || t == "NOJC" || t == "NOJD" || t == "NOJE" || t == "NOJF" ||
            t == "NOJG" || t == "NOJH" || t == "NOJI" || t == "NOJJ" || t == "NOJK" || t == "NOJL" ||
            t == "NOJM" || t == "NOJN" || t == "NOJO" || t == "NOJP" || t == "NOJQ" || t == "NOJR" ||
            t == "NOJS" || t == "NOJT" || t == "TRAV" || t == "NFRT" || t == "XOBL" || t == "NSHC") {
            for (int i = 0; i < 8; i++) (void)read_n3_value<float>();
            (void)read_n3_value<int32_t>();
            return true;
        }

        else if (t == "LKSN") {
            int32_t numSkinLists = read_n3_value<int32_t>();
            std::cout << "LKSN: numSkinLists=" << numSkinLists << std::endl;
            return true;
            }
        else if (t == "TXTS") {
            std::string key = read_n3_string();
            std::string val = read_n3_string();
            node->shader_textures[key] = val;
            return true;
        }
        else if (t == "TLFS") {
            std::string name = read_n3_string();
            float v = read_n3_value<float>();
            node->shader_parameters[name] = v;
            return true;
        }
        else if (t == "CEVS") {
            std::string name = read_n3_string();
            read_n3_value<float>();
            read_n3_value<float>();
            read_n3_value<float>();
            read_n3_value<float>();
            return true;

        }
        else if (t == "HSEM") {
            std::string mesh_res_id = read_n3_string();
            node->mesh_resource_id = mesh_res_id;
            std::cout << "MESH: " << mesh_res_id << std::endl;
            return true;
            }
        else if (t == "TNIS" || t == "SINT") { (void)read_n3_string(); (void)read_n3_value<int32_t>(); return true; }
        else if (t == "HSAC") {
            (void)read_n3_value<uint8_t>();
            return true;
        }

        else if (t == "EDFR") { return true; }
        else if (t == "RDHS" || t == "PTNM" || t == "SEMS") {
            (void)read_n3_string();
            return true;
        }
        else if (t == "IRGP") {
            (void)read_n3_value<uint32_t>();
            return true;
        }
        else if (t == "FKSN") {
            (void)read_n3_value<uint32_t>();
            return true;
        }
        else if (t == "GRFS") {
            uint32_t a0 = read_n3_value<uint32_t>();
            uint32_t n = read_n3_value<uint32_t>();
            if (n > 0 && n < 10000) {
                for (uint32_t i = 0; i < n; ++i) (void)read_n3_value<uint32_t>();
            }
            else {
                while (n3file.good()) {
                    std::string next = peek_fourcc(n3file, byteorder);
                    if (looks_like_tag(next)) break;
                    (void)read_n3_value<uint32_t>();
                }
            }
            return true;
        }

        return false;
    }



    bool parse_tag_transform(const std::string& t, std::shared_ptr<Node>& node) {
        if (t == "POSI") {
            float pos[4];
            for (int i = 0; i < 4; i++) pos[i] = read_n3_value<float>();
            node->position = { pos[0], pos[1], pos[2] };
            std::cout << " " << node->node_name << " POSI: " << pos[0] << ", " << pos[1] << ", " << pos[2] << std::endl;
            return true;
        }
        else if (t == "SCAL") {
            float scl[4];
            for (int i = 0; i < 4; i++) scl[i] = read_n3_value<float>();
            node->scale = { scl[0], scl[1], scl[2] };
            std::cout << " " << node->node_name << " SCAL: " << scl[0] << ", " << scl[1] << ", " << scl[2] << std::endl;
            return true;
        }
        else if (t == "ROTN") {
            float rot[4];
            for (int i = 0; i < 4; i++) rot[i] = read_n3_value<float>();
            node->rotation = { rot[0], rot[1], rot[2], rot[3] };
            std::cout << " " << node->node_name << " ROTN: " << rot[0] << ", " << rot[1] << ", " << rot[2] << ", " << rot[3] << std::endl;
            return true;
        }
        return false;
    }
    bool parse_tag_state(const std::string& t, std::shared_ptr<Node>& node) {
        if (t == "MNMT") { read_n3_string(); return true; }
        else if (t == "MATE") { read_n3_string(); return true; }
        else if (t == "STXT") {
            std::string k = read_n3_string();
            std::string v = read_n3_string();
            node->shader_textures[k] = v;
            currentSlot = infer_slot_from_key(k);
            if (!node->uvXformBySlot.count(currentSlot)) node->uvXformBySlot[currentSlot] = std::array<float, 4>{ 1.0f,1.0f,0.0f,0.0f };
            if (!node->uvSetBySlot.count(currentSlot)) node->uvSetBySlot[currentSlot] = 0;
            return true;
        }
        return false;
    }

    
public:
    bool parse_file(const std::string& path) {
        n3file.open(path, std::ios::binary);
        if (!n3file.is_open()) {
            report("ERROR", "Could not open file: " + path);
            return false;
        }

        char hdr[4];
        n3file.read(hdr, 4);
        std::string h(hdr, 4);

        if (h == "NEB3") byteorder = "little";
        else if (h == "3BEN") byteorder = "big";

        uint32_t ver = read_n3_value<uint32_t>();
        std::cout << "n3 Version: " << ver << std::endl;
        n3version = ver;

        bool done = false;
        std::shared_ptr<Node> current = nullptr;
        std::vector<std::shared_ptr<Node>> stack;

        while (!done && n3file.good()) {
            std::streampos pos_before = n3file.tellg();
            std::string tag;
            try {
                tag = read_n3_fourcc();
              //  std::cout << "Position: 0x" << std::hex << pos_before << std::dec
              //      << ", Tag: '" << tag << "'" << std::endl;
            }
            catch (std::runtime_error& e) {
                std::cout << "Exception at position 0x" << std::hex << pos_before << std::dec
                    << ": " << e.what() << std::endl;
                if (resync_to_valid()) continue;
                break;
            }

            if (ValidFourCC().count(tag) == 0) {
                errout << "Non-FourCC or garbage '" << tag << "' at 0x"
                    << std::hex << pos_before << std::dec << " -> trying to skip/resync\n";

                n3file.clear();
                n3file.seekg(pos_before);
                if (generic_skip()) continue;
                if (resync_to_valid()) continue;
                report("ERROR", "Cannot resync at position 0x" +
                    std::to_string(static_cast<long long>(pos_before)));
                break;
            }

            if (tag == "LDM>") {
                n3modeltype = read_n3_fourcc();
                n3modelname = read_n3_string();
                std::cout << "LDM> fourcc: '" << n3modeltype << "'" << std::endl;
                std::cout << "LDM> name: '" << n3modelname << "'" << std::endl;

                auto n = std::make_shared<Node>();
                n->node_type = n3modeltype;
                n->node_name = n3modelname;
                n->node_parent = nullptr;
                n3node_list.push_back(n);
                current = n;
                stack.push_back(n);
                currentSlot = 0;
            }
            else if (tag == "<LDM") {
                done = true;
                current = nullptr;
            }
            else if (tag == "DNM>") {
                std::string node4 = read_n3_fourcc();
                std::string nodeName = read_n3_string();

                auto n = std::make_shared<Node>();
                n->node_name = nodeName;
                n->node_type = node4;
                n->node_parent = current;
                std::cout << "DNM> " << n->node_type << " - " << n->node_name << std::endl;
                if (current) current->node_children.push_back(n);
                n3node_list.push_back(n);
                current = n;
                stack.push_back(n);
                currentSlot = 0;
            }
            else if (tag == "<DNM") {
                if (!stack.empty()) stack.pop_back();
                current = stack.empty() ? nullptr : stack.back();
            }
            else if (tag == "EOF_" || tag == "_FOE") {
                done = true;
                break;
            }
            else {
                if (!current) {
                    report("ERROR", "Tag without active node: " + tag);
                    return false;
                }
                if (!parse_tag_shape(tag, current) &&
                    !parse_tag_state(tag, current) &&
                    !parse_tag_transform(tag, current) &&
                    !skip4cc(tag, current)) {
                    if (!generic_skip()) {
                        report("ERROR", "Unknown tag '" + tag + "' at position 0x" +
                            std::to_string(static_cast<long long>(pos_before)));
                        return false;
                    }
                }
            }

            std::streampos pos_after = n3file.tellg();
           // std::cout << "After processing, position: 0x" << std::hex << pos_after << std::dec << std::endl;
        }

        n3file.close();
        return true;
    }
};


void computeHierarchyMatrices(const std::vector<std::shared_ptr<Node>>& nodes) {
    std::function<void(const std::shared_ptr<Node>&)> dfs = [&](const std::shared_ptr<Node>& n) {
        if (n->node_parent) {
            float tmp[16];
            multiplyMatrices(tmp, n->node_parent->worldMatrix, n->localMatrix);
            memcpy(n->worldMatrix, tmp, sizeof(tmp));
        }
        else {
            memcpy(n->worldMatrix, n->localMatrix, sizeof(n->localMatrix));
        }
        for (auto& c : n->node_children) dfs(c);
        };
    for (auto& n : nodes) if (!n->node_parent) dfs(n);
}


GLuint LoadNebTexSmart(const std::string& path) {
    if (path.find("system/white.dds") != std::string::npos) return gWhiteTex;
    if (path.find("system/black.dds") != std::string::npos) return gBlackTex;
    if (path.find("system/nobump.dds") != std::string::npos) return gFlatNormalTex;
    return LoadDDS(path);
}


static GLuint loadNebTex(const std::unordered_map<std::string, std::string>& m, const char* key) {
    auto it = m.find(key);
    if (it == m.end()) return 0;
    std::string t = it->second;
    if (t.rfind("tex:", 0) == 0) t = t.substr(4);
    return LoadNebTexSmart(std::string("C:/drasa_online/work/textures/") + t + ".dds");
}

void bindNodeTextures(const Node& node, GLint uDiff, GLint uSpec, GLint uNorm, GLint uEmis) {
    GLuint td = loadNebTex(node.shader_textures, "DiffMap0");
    GLuint ts = loadNebTex(node.shader_textures, "SpecMap0");
    GLuint tn = loadNebTex(node.shader_textures, "BumpMap0");
    GLuint te = loadNebTex(node.shader_textures, "EmissiveMap0");

    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, td); glUniform1i(uDiff, 0);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, ts); glUniform1i(uSpec, 1);
    glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, tn); glUniform1i(uNorm, 2);
    glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, te); glUniform1i(uEmis, 3);
}



void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse) {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    double dx = xpos - lastX;
    double dy = ypos - lastY;
    lastX = xpos;
    lastY = ypos;

    if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        yaw += (float)dx * 0.01f;
        pitch -= (float)dy * 0.01f;
        pitch = clamp(pitch, -1.5f, 1.5f);
    }
    else {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }
}

void InitGrid(float gridSize, float gridStep) {
    struct Vertex { float x, y, z; };
    std::vector<Vertex> gridVertices;

    for (float i = -gridSize; i <= gridSize; i += gridStep) {
        gridVertices.push_back({ i, 0.0f, -gridSize });
        gridVertices.push_back({ i, 0.0f,  gridSize });
        gridVertices.push_back({ -gridSize, 0.0f, i });
        gridVertices.push_back({ gridSize, 0.0f, i });
    }
    gridVertexCount = (GLsizei)gridVertices.size();

    glGenVertexArrays(1, &gridVAO);
    glGenBuffers(1, &gridVBO);

    glBindVertexArray(gridVAO);
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
    glBufferData(GL_ARRAY_BUFFER, gridVertexCount * sizeof(Vertex), gridVertices.data(), GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);

    glBindVertexArray(0);
}

void DrawGrid(GLuint shaderProgram, const float* view, const float* proj) {
    glUseProgram(shaderProgram);

    float model[16] = {
        1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1
    };

    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, model);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, view);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, proj);
    glUniform3f(glGetUniformLocation(shaderProgram, "objectColor"), 0.15f, 0.15f, 0.15f);

    glBindVertexArray(gridVAO);
    glDrawArrays(GL_LINES, 0, gridVertexCount);
    glBindVertexArray(0);

    glUseProgram(0);
}


void CleanupGrid() {
    glDeleteVertexArrays(1, &gridVAO);
    glDeleteBuffers(1, &gridVBO);
}

void processInput(GLFWwindow* window)
{
    float speed = 0.0001f * distance;
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        speed *= 1.0f;

    float rightX = cos(yaw);
    float rightZ = sin(yaw);

    float forwardX = -sin(yaw) * cos(pitch);
    float forwardY = -sin(pitch);
    float forwardZ = -cos(yaw) * cos(pitch);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        targetX += forwardX * speed;
        targetY += forwardY * speed;
        targetZ += forwardZ * speed;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        targetX -= forwardX * speed;
        targetY -= forwardY * speed;
        targetZ -= forwardZ * speed;
    }

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        targetX -= rightX * speed;
        targetZ -= rightZ * speed;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        targetX += rightX * speed;
        targetZ += rightZ * speed;
    }

    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) targetY -= speed;
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) targetY += speed;

    int cState = glfwGetKey(window, GLFW_KEY_C);
    if (cState == GLFW_PRESS && !gCPrev) {
        gShowDDS = !gShowDDS;
    }
    gCPrev = (cState == GLFW_PRESS);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    distance -= static_cast<float>(yoffset) * 0.5f;
    if (distance < 1.0f) distance = 1.0f;
    if (distance > 50.0f) distance = 50.0f;
}

float gridVertices[] = {
    -5.0f, 0.0f, -5.0f,
     5.0f, 0.0f, -5.0f,
     5.0f, 0.0f,  5.0f,
    -5.0f, 0.0f,  5.0f,
    -5.0f, 0.0f, -5.0f,
     5.0f, 0.0f, -5.0f,
     5.0f, 0.0f,  5.0f,
    -5.0f, 0.0f,  5.0f
};

enum VertexComponent {
    Coord = 1 << 0,
    Normal = 1 << 1,
    NormalUB4N = 1 << 2,
    Uv0 = 1 << 3,
    Uv0S2 = 1 << 4,
    Uv1 = 1 << 5,
    Uv1S2 = 1 << 6,
    Uv2 = 1 << 7,
    Uv2S2 = 1 << 8,
    Uv3 = 1 << 9,
    Uv3S2 = 1 << 10,
    Color = 1 << 11,
    ColorUB4N = 1 << 12,
    Tangent = 1 << 13,
    TangentUB4N = 1 << 14,
    Binormal = 1 << 15,
    BinormalUB4N = 1 << 16,
    Weights = 1 << 17,
    WeightsUB4N = 1 << 18,
    JIndices = 1 << 19,
    JIndicesUB4 = 1 << 20
};

static int GetComponentSize(uint32_t comp) {
    switch (comp) {
    case Coord: return 12;
    case Normal: return 12;
    case NormalUB4N: return 4;
    case Uv0: return 8;
    case Uv0S2: return 4;
    case Uv1: return 8;
    case Uv1S2: return 4;
    case Uv2: return 8;
    case Uv2S2: return 4;
    case Uv3: return 8;
    case Uv3S2: return 4;
    case Color: return 16;
    case ColorUB4N: return 4;
    case Tangent: return 12;
    case TangentUB4N: return 4;
    case Binormal: return 12;
    case BinormalUB4N: return 4;
    case Weights: return 16;
    case WeightsUB4N: return 4;
    case JIndices: return 16;
    case JIndicesUB4: return 4;
    default: return 0;
    }
}

/*
static void RecalcNormals(Mesh& out) {
    for (auto& v : out.verts) { v.nx = v.ny = v.nz = 0.0f; }
    for (size_t i = 0; i + 2 < out.idx.size(); i += 3) {
        ObjVertex& v0 = out.verts[out.idx[i + 0]];
        ObjVertex& v1 = out.verts[out.idx[i + 1]];
        ObjVertex& v2 = out.verts[out.idx[i + 2]];

        float ux = v1.px - v0.px, uy = v1.py - v0.py, uz = v1.pz - v0.pz;
        float vx = v2.px - v0.px, vy = v2.py - v0.py, vz = v2.pz - v0.pz;
        float nx = uy * vz - uz * vy;
        float ny = uz * vx - ux * vz;
        float nz = ux * vy - uy * vx;
        float len = sqrtf(nx * nx + ny * ny + nz * nz);
        if (len > 1e-6f) { nx /= len; ny /= len; nz /= len; }

        v0.nx += nx; v0.ny += ny; v0.nz += nz;
        v1.nx += nx; v1.ny += ny; v1.nz += nz;
        v2.nx += nx; v2.ny += ny; v2.nz += nz;
    }
    for (auto& v : out.verts) {
        float len = sqrtf(v.nx * v.nx + v.ny * v.ny + v.nz * v.nz);
        if (len > 1e-6f) { v.nx /= len; v.ny /= len; v.nz /= len; }
        else { v.nx = 0; v.ny = 0; v.nz = 1; }
    }
}*/


static inline float ub_to_n11(uint8_t b) { return (float(b) * (1.0f / 255.0f)) * 2.0f - 1.0f; } static inline float ub_to_u01(uint8_t b) { return float(b) * (1.0f / 255.0f); } static inline float s_to_n11(int16_t s) { return s < 0 ? float(s) / 32768.0f : float(s) / 32767.0f; } static inline float s_to_u01(int16_t s) { return (float(s) + 32768.0f) / 65535.0f; }
static inline float s2_to_uv(int16_t s) { return float(s) / 8191.0f; }

static inline uint32_t bit(uint32_t x) { return 1u << x; }

static void SetupVAO(Mesh& out) {
    if (out.vao == 0) glGenVertexArrays(1, &out.vao);
    if (out.vbo == 0) glGenBuffers(1, &out.vbo);
    if (out.ebo == 0) glGenBuffers(1, &out.ebo);
    glBindVertexArray(out.vao); glBindBuffer(GL_ARRAY_BUFFER, out.vbo);
    glBufferData(GL_ARRAY_BUFFER, out.verts.size() * sizeof(ObjVertex), out.verts.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, out.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, out.idx.size() * sizeof(uint32_t), out.idx.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, px));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, nx));
    glEnableVertexAttribArray(2); glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, u0));
    glEnableVertexAttribArray(3); glVertexAttribPointer(3, 4, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, w0));
    glEnableVertexAttribArray(4); glVertexAttribIPointer(4, 4, GL_UNSIGNED_BYTE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, j0));
   
    glEnableVertexAttribArray(5);

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, tx));

    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, bx));

    glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, u1));

    glBindVertexArray(0);
}
static inline float fps2float(int16_t s) {
    return static_cast<float>(s) / 8191.0f;
}

/*
bool LoadNVX2_AndNormalizeUnit(const std::string& path, Mesh& out) {
    std::cout << "[NVX2] open: " << path << "\n";
    std::ifstream f(path, std::ios::binary);
    if (!f) { std::cerr << "[NVX2] cannot open\n"; return false; }

    char magicChars[4] = { 0,0,0,0 };
    f.read(magicChars, 4);
    if (!f) { std::cerr << "[NVX2] failed reading magic\n"; return false; }
    std::string magicStr(magicChars, 4);
    std::cout << "[NVX2] magic: '" << magicStr << "' ("
        << std::hex << std::showbase
        << (int)(unsigned char)magicChars[0] << " "
        << (int)(unsigned char)magicChars[1] << " "
        << (int)(unsigned char)magicChars[2] << " "
        << (int)(unsigned char)magicChars[3] << std::dec << ")\n";

    uint32_t numGroups = 0, numVertices = 0, vertexWidthFloats = 0;
    uint32_t numTrianglesOrIndices = 0, numEdges = 0, compMask = 0;

    f.read((char*)&numGroups, 4);
    f.read((char*)&numVertices, 4);
    f.read((char*)&vertexWidthFloats, 4);
    f.read((char*)&numTrianglesOrIndices, 4);
    f.read((char*)&numEdges, 4);
    f.read((char*)&compMask, 4);

    std::cout << "[NVX2] header numGroups=" << numGroups
        << " numVertices=" << numVertices
        << " vertexWidthFloats=" << vertexWidthFloats
        << " triOrIdx=" << numTrianglesOrIndices
        << " numEdges=" << numEdges
        << " compMask=0x" << std::hex << compMask << std::dec << "\n";
    std::cout << "[NVX2] tellg after header: " << (long long)f.tellg() << "\n";

    out.groups.clear();
    out.groups.resize(numGroups);
    if (numGroups > 0) {
        f.read((char*)out.groups.data(), numGroups * sizeof(Nvx2Group));
        if (!f) { std::cerr << "[NVX2] failed reading groups\n"; return false; }
    }
    std::cout << "[NVX2] read groups: " << numGroups << " tellg=" << (long long)f.tellg() << "\n";
    for (uint32_t gi = 0; gi < numGroups; ++gi) {
        const auto& g = out.groups[gi];
        std::cout << "[NVX2] group[" << gi << "] v=(" << g.firstVertex << "," << g.numVertices
            << ") tri=(" << g.firstTriangle << "," << g.numTriangles
            << ") edge=(" << g.firstEdge << "," << g.numEdges << ")\n";
    }

    const uint32_t strideBytesHeader = vertexWidthFloats * 4;
    std::cout << "[NVX2] strideBytesHeader=" << strideBytesHeader << "\n";

    const uint32_t order[] = {
        Coord, Normal, NormalUB4N,
        Uv0, Uv0S2, Uv1, Uv1S2, Uv2, Uv2S2, Uv3, Uv3S2,
        Color, ColorUB4N,
        Tangent, TangentUB4N,
        Binormal, BinormalUB4N,
        Weights, WeightsUB4N,
        JIndices, JIndicesUB4
    };

    std::unordered_map<uint32_t, int> offsets;
    int currentOffset = 0;
    auto addIf = [&](uint32_t comp, const char* name) {
        if (compMask & comp) {
            int sz = GetComponentSize(comp);
            offsets[comp] = currentOffset;
            std::cout << "[NVX2] enable " << name << " off=" << currentOffset << " size=" << sz << "\n";
            currentOffset += sz;
        }
        else {
            std::cout << "[NVX2] disable " << name << "\n";
        }
        };

    addIf(Coord, "Coord");
    addIf(Normal, "Normal");
    addIf(NormalUB4N, "NormalUB4N");
    addIf(Uv0, "Uv0");
    addIf(Uv0S2, "Uv0S2");
    addIf(Uv1, "Uv1");
    addIf(Uv1S2, "Uv1S2");
    addIf(Uv2, "Uv2");
    addIf(Uv2S2, "Uv2S2");
    addIf(Uv3, "Uv3");
    addIf(Uv3S2, "Uv3S2");
    addIf(Color, "Color");
    addIf(ColorUB4N, "ColorUB4N");
    addIf(Tangent, "Tangent");
    addIf(TangentUB4N, "TangentUB4N");
    addIf(Binormal, "Binormal");
    addIf(BinormalUB4N, "BinormalUB4N");
    addIf(Weights, "Weights");
    addIf(WeightsUB4N, "WeightsUB4N");
    addIf(JIndices, "JIndices");
    addIf(JIndicesUB4, "JIndicesUB4");

    std::cout << "[NVX2] computedStride=" << currentOffset << "\n";
    if ((uint32_t)currentOffset != strideBytesHeader) {
        std::cerr << "[NVX2] stride mismatch header=" << strideBytesHeader << " computed=" << currentOffset << "\n";
    }

    std::vector<uint8_t> vbuf;
    vbuf.resize(size_t(numVertices) * size_t(strideBytesHeader));
    if (!vbuf.empty()) {
        f.read((char*)vbuf.data(), vbuf.size());
        if (!f) { std::cerr << "[NVX2] failed reading vertex buffer bytes=" << vbuf.size() << "\n"; return false; }
    }
    std::cout << "[NVX2] read VBUF bytes=" << vbuf.size() << " tellg=" << (long long)f.tellg() << "\n";

    out.verts.clear();
    out.verts.resize(numVertices);

    bool hasUv0 = offsets.count(Uv0) || offsets.count(Uv0S2);
    bool hasUv1 = offsets.count(Uv1) || offsets.count(Uv1S2);
    std::cout << "[NVX2] hasUv0=" << (hasUv0 ? 1 : 0) << " hasUv1=" << (hasUv1 ? 1 : 0) << "\n";

    for (uint32_t i = 0; i < numVertices; i++) {
        const uint8_t* base = vbuf.data() + size_t(i) * size_t(strideBytesHeader);
        ObjVertex v{};
        v.nz = 1.0f; v.cr = v.cg = v.cb = v.ca = 1.0f; v.w0 = 1.0f;

        if (offsets.count(Coord)) {
            const float* p = (const float*)(base + offsets[Coord]);
            v.px = p[0]; v.py = p[1]; v.pz = p[2];
        }
        if (offsets.count(Normal)) {
            const float* n = (const float*)(base + offsets[Normal]);
            v.nx = n[0]; v.ny = n[1]; v.nz = n[2];
        }

        if (offsets.count(Tangent)) {
            const float* t = (const float*)(base + offsets[Tangent]);
            v.tx = t[0]; v.ty = t[1]; v.tz = t[2];
        }
        else if (offsets.count(TangentUB4N)) {
            const uint8_t* t = base + offsets[TangentUB4N];
            v.tx = ub_to_n11(t[0]); v.ty = ub_to_n11(t[1]); v.tz = ub_to_n11(t[2]);
        }

        if (offsets.count(Binormal)) {
            const float* b = (const float*)(base + offsets[Binormal]);
            v.bx = b[0]; v.by = b[1]; v.bz = b[2];
        }
        else if (offsets.count(BinormalUB4N)) {
            const uint8_t* b = base + offsets[BinormalUB4N];
            v.bx = ub_to_n11(b[0]); v.by = ub_to_n11(b[1]); v.bz = ub_to_n11(b[2]);
        }

        if (offsets.count(Uv0)) {
            const float* u = (const float*)(base + offsets[Uv0]);
            v.u0 = u[0]; v.v0 = u[1];
        }
        if (offsets.count(Uv0S2)) {
            const int16_t* u = (const int16_t*)(base + offsets[Uv0S2]);
            v.u0 = fps2float(u[0]);
            v.v0 = fps2float(u[1]);
        }

        if (offsets.count(Uv1)) {
            const float* u = (const float*)(base + offsets[Uv1]);
            v.u1 = u[0]; v.v1 = 1.0f - u[1];
        }
        if (offsets.count(Uv1S2)) {
            const int16_t* u = (const int16_t*)(base + offsets[Uv1S2]);
            v.u1 = fps2float(u[0]);
            v.v1 = fps2float(u[1]);
        }


        if (offsets.count(Color)) {
            const float* c = (const float*)(base + offsets[Color]);
            v.cr = c[0]; v.cg = c[1]; v.cb = c[2]; v.ca = c[3];
        }
        else if (offsets.count(ColorUB4N)) {
            const uint8_t* c = base + offsets[ColorUB4N];
            v.cr = ub_to_u01(c[0]); v.cg = ub_to_u01(c[1]); v.cb = ub_to_u01(c[2]); v.ca = ub_to_u01(c[3]);
        }

        if (offsets.count(Weights)) {
            const float* w = (const float*)(base + offsets[Weights]);
            v.w0 = w[0]; v.w1 = w[1]; v.w2 = w[2]; v.w3 = w[3];
        }
        else if (offsets.count(WeightsUB4N)) {
            const uint8_t* w = (base + offsets[WeightsUB4N]);
            v.w0 = ub_to_u01(w[0]); v.w1 = ub_to_u01(w[1]); v.w2 = ub_to_u01(w[2]); v.w3 = ub_to_u01(w[3]);
            float s = v.w0 + v.w1 + v.w2 + v.w3; if (s > 0) { float inv = 1.0f / s; v.w0 *= inv; v.w1 *= inv; v.w2 *= inv; v.w3 *= inv; }
        }

        if (offsets.count(JIndices)) {
            const float* jf = (const float*)(base + offsets[JIndices]);
            v.j0 = (uint8_t)jf[0]; v.j1 = (uint8_t)jf[1]; v.j2 = (uint8_t)jf[2]; v.j3 = (uint8_t)jf[3];
        }
        else if (offsets.count(JIndicesUB4)) {
            const uint8_t* jb = (base + offsets[JIndicesUB4]);
            v.j0 = jb[0]; v.j1 = jb[1]; v.j2 = jb[2]; v.j3 = jb[3];
        }

        out.verts[i] = v;

        if (i < 3) {
            std::cout << "[NVX2] v[" << i << "] P(" << v.px << "," << v.py << "," << v.pz << ") N("
                << v.nx << "," << v.ny << "," << v.nz << ") UV0(" << v.u0 << "," << v.v0
                << ") UV1(" << v.u1 << "," << v.v1 << ") W(" << v.w0 << "," << v.w1 << ","
                << v.w2 << "," << v.w3 << ") J(" << (int)v.j0 << "," << (int)v.j1 << ","
                << (int)v.j2 << "," << (int)v.j3 << ")\n";
        }
    }
    std::cout << "[NVX2] verts=" << out.verts.size() << "\n";

    uint32_t totalTris = 0;
    for (auto& g : out.groups) totalTris += g.numTriangles;
    uint32_t expectedIndexCount = totalTris ? totalTris * 3 : numTrianglesOrIndices;
    std::cout << "[NVX2] totalTris=" << totalTris << " expectedIndexCount=" << expectedIndexCount << "\n";
    std::cout << "[NVX2] tellg before IDX: " << (long long)f.tellg() << "\n";

    out.idx.clear();
    if (expectedIndexCount > 0) {
        std::vector<uint16_t> tmp(expectedIndexCount);
        f.read((char*)tmp.data(), expectedIndexCount * sizeof(uint16_t));
        if (!f) {
            std::cerr << "[NVX2] failed reading indices count=" << expectedIndexCount << "\n";
            return false;
        }
        out.idx.resize(expectedIndexCount);
        for (uint32_t i = 0; i < expectedIndexCount; i++) out.idx[i] = tmp[i];
        std::cout << "[NVX2] read indices u16 count=" << out.idx.size() << " tellg=" << (long long)f.tellg() << "\n";
    }
    else {
        std::cout << "[NVX2] no indices\n";
    }

    if (!out.idx.empty()) {
        std::cout << "[NVX2] idx[0.." << std::min<size_t>(5, out.idx.size()) - 1 << "]:";
        for (size_t i = 0; i < std::min<size_t>(5, out.idx.size()); ++i) std::cout << " " << out.idx[i];
        std::cout << "\n";
    }
    RecalcNormals(out);

    SetupVAO(out);
    std::cout << "[NVX2] VAO=" << out.vao << " VBO=" << out.vbo << " EBO=" << out.ebo
        << " verts=" << out.verts.size() << " idx=" << out.idx.size() << "\n";

    return true;
}
*/

bool LoadNVX2(const std::string& path, Mesh& out) {
    const float px = 0.0f, py = 0.0f, pz = 0.0f;
    const float rx = 0.0f, ry = 0.0f, rz = 0.0f;
    const float sx = 1.0f, sy = 1.0f, sz = 1.0f;

    std::ifstream f(path, std::ios::binary);
    if (!f) { errout << "[NVX2] cannot open " << path << " "; return false; }

    char magicChars[4] = { 0,0,0,0 };
    f.read(magicChars, 4);
    if (!f) { errout << "[NVX2] failed reading magic\n"; return false; }
    std::string magicStr(magicChars, 4);
    std::cout << "[NVX2] magic: '" << magicStr << "' ("
        << std::hex << std::showbase
        << (int)(unsigned char)magicChars[0] << " "
        << (int)(unsigned char)magicChars[1] << " "
        << (int)(unsigned char)magicChars[2] << " "
        << (int)(unsigned char)magicChars[3] << std::dec << ")\n";

    uint32_t numGroups, numVertices, vertexWidthFloats, numTrianglesOrIndices, numEdges, compMask;
    f.read((char*)&numGroups, 4);
    f.read((char*)&numVertices, 4);
    f.read((char*)&vertexWidthFloats, 4);
    f.read((char*)&numTrianglesOrIndices, 4);
    f.read((char*)&numEdges, 4);
    f.read((char*)&compMask, 4);

    out.groups.clear();
    out.groups.resize(numGroups);
    if (numGroups > 0) {
        f.read((char*)out.groups.data(), numGroups * sizeof(Nvx2Group));
        if (!f) return false;
    }

    const uint32_t strideBytesHeader = vertexWidthFloats * 4;

    std::unordered_map<uint32_t, int> offsets;
    int currentOffset = 0;
    auto addIf = [&](uint32_t comp) {
        if (compMask & comp) {
            offsets[comp] = currentOffset;
            currentOffset += GetComponentSize(comp);
        }
        };

    addIf(Coord); addIf(Normal); addIf(NormalUB4N);
    addIf(Uv0); addIf(Uv0S2); addIf(Uv1); addIf(Uv1S2);
    addIf(Uv2); addIf(Uv2S2); addIf(Uv3); addIf(Uv3S2);
    addIf(Color); addIf(ColorUB4N);
    addIf(Tangent); addIf(TangentUB4N);
    addIf(Binormal); addIf(BinormalUB4N);
    addIf(Weights); addIf(WeightsUB4N);
    addIf(JIndices); addIf(JIndicesUB4);

    std::vector<uint8_t> vbuf(size_t(numVertices) * size_t(strideBytesHeader));
    if (!vbuf.empty()) {
        f.read((char*)vbuf.data(), vbuf.size());
        if (!f) return false;
    }

    out.verts.clear();
    out.verts.resize(numVertices);

    float cx = cosf(rx), sx_r = sinf(rx);
    float cy = cosf(ry), sy_r = sinf(ry);
    float cz = cosf(rz), sz_r = sinf(rz);

    for (uint32_t i = 0; i < numVertices; i++) {
        const uint8_t* base = vbuf.data() + size_t(i) * size_t(strideBytesHeader);
        ObjVertex v{};
        v.nz = 1.0f; v.cr = v.cg = v.cb = v.ca = 1.0f; v.w0 = 1.0f;

        if (offsets.count(Coord)) {
            const float* p = (const float*)(base + offsets[Coord]);
            float x = p[0] * sx, y = p[1] * sy, z = p[2] * sz;

            float x1 = x, y1 = y * cy - z * sy_r, z1 = y * sy_r + z * cy;
            float x2 = x1 * cx + z1 * sx_r, y2 = y1, z2 = -x1 * sx_r + z1 * cx;
            v.px = x2 * cz - y2 * sz_r + px;
            v.py = x2 * sz_r + y2 * cz + py;
            v.pz = z2 + pz;
        }

        if (offsets.count(Normal)) {
            const float* n = (const float*)(base + offsets[Normal]);
            v.nx = n[0]; v.ny = n[1]; v.nz = n[2];
        }

        if (offsets.count(Tangent)) {
            const float* t = (const float*)(base + offsets[Tangent]);
            v.tx = t[0]; v.ty = t[1]; v.tz = t[2];
        }
        else if (offsets.count(TangentUB4N)) {
            const uint8_t* t = base + offsets[TangentUB4N];
            v.tx = ub_to_n11(t[0]); v.ty = ub_to_n11(t[1]); v.tz = ub_to_n11(t[2]);
        }

        if (offsets.count(Binormal)) {
            const float* b = (const float*)(base + offsets[Binormal]);
            v.bx = b[0]; v.by = b[1]; v.bz = b[2];
        }
        else if (offsets.count(BinormalUB4N)) {
            const uint8_t* b = base + offsets[BinormalUB4N];
            v.bx = ub_to_n11(b[0]); v.by = ub_to_n11(b[1]); v.bz = ub_to_n11(b[2]);
        }

        if (offsets.count(Uv0)) {
            const float* u = (const float*)(base + offsets[Uv0]);
            v.u0 = u[0]; v.v0 = u[1];
        }
        if (offsets.count(Uv0S2)) {
            const int16_t* u = (const int16_t*)(base + offsets[Uv0S2]);
            v.u0 = fps2float(u[0]); v.v0 = fps2float(u[1]);
        }

        if (offsets.count(Uv1)) {
            const float* u = (const float*)(base + offsets[Uv1]);
            v.u1 = u[0]; v.v1 = 1.0f - u[1];
        }
        if (offsets.count(Uv1S2)) {
            const int16_t* u = (const int16_t*)(base + offsets[Uv1S2]);
            v.u1 = fps2float(u[0]); v.v1 = fps2float(u[1]);
        }

        if (offsets.count(Color)) {
            const float* c = (const float*)(base + offsets[Color]);
            v.cr = c[0]; v.cg = c[1]; v.cb = c[2]; v.ca = c[3];
        }
        else if (offsets.count(ColorUB4N)) {
            const uint8_t* c = base + offsets[ColorUB4N];
            v.cr = ub_to_u01(c[0]); v.cg = ub_to_u01(c[1]); v.cb = ub_to_u01(c[2]); v.ca = ub_to_u01(c[3]);
        }

        if (offsets.count(Weights)) {
            const float* w = (const float*)(base + offsets[Weights]);
            v.w0 = w[0]; v.w1 = w[1]; v.w2 = w[2]; v.w3 = w[3];
        }
        else if (offsets.count(WeightsUB4N)) {
            const uint8_t* w = (base + offsets[WeightsUB4N]);
            v.w0 = ub_to_u01(w[0]); v.w1 = ub_to_u01(w[1]); v.w2 = ub_to_u01(w[2]); v.w3 = ub_to_u01(w[3]);
            float s = v.w0 + v.w1 + v.w2 + v.w3;
            if (s > 0) { float inv = 1.0f / s; v.w0 *= inv; v.w1 *= inv; v.w2 *= inv; v.w3 *= inv; }
        }

        if (offsets.count(JIndices)) {
            const float* jf = (const float*)(base + offsets[JIndices]);
            v.j0 = (uint8_t)jf[0]; v.j1 = (uint8_t)jf[1]; v.j2 = (uint8_t)jf[2]; v.j3 = (uint8_t)jf[3];
        }
        else if (offsets.count(JIndicesUB4)) {
            const uint8_t* jb = (base + offsets[JIndicesUB4]);
            v.j0 = jb[0]; v.j1 = jb[1]; v.j2 = jb[2]; v.j3 = jb[3];
        }

        out.verts[i] = v;
    }

    uint32_t totalTris = 0;
    for (auto& g : out.groups) totalTris += g.numTriangles;
    uint32_t expectedIndexCount = totalTris ? totalTris * 3 : numTrianglesOrIndices;

    out.idx.clear();
    if (expectedIndexCount > 0) {
        std::vector<uint16_t> tmp(expectedIndexCount);
        f.read((char*)tmp.data(), expectedIndexCount * sizeof(uint16_t));
        if (!f) return false;
        out.idx.resize(expectedIndexCount);
        for (uint32_t i = 0; i < expectedIndexCount; i++) out.idx[i] = tmp[i];
    }

    SetupVAO(out);
    return true;
}

const char* vertexShaderSource = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 aNormal;
layout(location=2) in vec2 aUV0; 
layout(location=5) in vec2 aUV1;

uniform mat4 model, view, projection;

out vec3 vPos;
out vec3 vNormal;
out vec2 vUV0;
out vec2 vUV1;

void main() {
    vec4 wp = model * vec4(aPos, 1.0);
    vPos = wp.xyz;
    vNormal = mat3(transpose(inverse(model))) * aNormal;
    vUV0 = aUV0;
    vUV1 = aUV1;
    gl_Position = projection * view * wp;
}

)";

const char* fragmentShaderSource = R"(
#version 330 core
in vec3 vPos;
in vec3 vNormal;
in vec2 vUV0;
in vec2 vUV1;

out vec4 FragColor;

uniform sampler2D DiffMap;
uniform sampler2D SpecMap;
uniform sampler2D BumpMap;
uniform sampler2D EmissiveMap;

uniform vec3 CameraPos;

uniform vec3 SunDir;
uniform vec3 SunColor;
uniform float SunIntensity;
uniform vec3 AmbientColor;

uniform float SpecularPower;
uniform float SpecularIntensity;
uniform float EmissiveIntensity;

uniform vec4 UvXform[4]; 
uniform int  UvSet[4]; 
uniform int  HasBump;
uniform int  HasSpec;

vec2 pickUV(int slot) {
    vec2 base = (UvSet[slot] == 1) ? vUV1 : vUV0;
    return base * UvXform[slot].xy + UvXform[slot].zw;
}

void main() {
    vec2 uvD = pickUV(0);
    vec4 diffTex = texture(DiffMap, uvD);
    if (diffTex.a < 0.5) discard;

    vec3 N = normalize(vNormal);
    if (HasBump == 1) {
        vec3 mapN = texture(BumpMap, pickUV(2)).xyz * 2.0 - 1.0;
        N = normalize(mix(N, mapN, 1.0));
    }

    vec3 V = normalize(CameraPos - vPos);
    vec3 L = normalize(-SunDir);
    vec3 H = normalize(L + V);

    float ndl = max(dot(N, L), 0.0);
    float ndh = max(dot(N, H), 0.0);

    vec3 albedo = diffTex.rgb;
    vec3 specC = (HasSpec == 1) ? texture(SpecMap, pickUV(1)).rgb : vec3(0.0);
    float spec = (ndl > 0.0) ? pow(ndh, SpecularPower) * SpecularIntensity : 0.0;
    vec3 emissive = texture(EmissiveMap, pickUV(3)).rgb * EmissiveIntensity;

    vec3 ambient = AmbientColor * albedo;
    vec3 diffuse = SunColor * SunIntensity * ndl * albedo;
    vec3 specular = SunColor * SunIntensity * spec * specC;

    vec3 color = ambient + diffuse + specular + emissive;

    FragColor = vec4(color, diffTex.a);
})";


static glm::mat4 composeTRS(const glm::vec3& p, const glm::vec3& eulerDeg, const glm::vec3& s) {
    glm::mat4 m(1.0f);
    m = glm::translate(m, p);
    m = glm::rotate(m, glm::radians(eulerDeg.z), glm::vec3(0, 0, 1));
    m = glm::rotate(m, glm::radians(eulerDeg.y), glm::vec3(0, 1, 0));
    m = glm::rotate(m, glm::radians(eulerDeg.x), glm::vec3(1, 0, 0));
    m = glm::scale(m, s);
    return m;
}

static void matCopy(float* dst16, const glm::mat4& m) {
    const float* p = glm::value_ptr(m);
    for (int i = 0; i < 16; ++i) dst16[i] = p[i];
}

static glm::mat4 matFromPtr(const float* src16) {
    glm::mat4 m(1.0f);
    float* p = glm::value_ptr(m);
    for (int i = 0; i < 16; ++i) p[i] = src16[i];
    return m;
}

static void setUvIdentity(float uv[4]) { uv[0] = 1.0f; uv[1] = 1.0f; uv[2] = 0.0f; uv[3] = 0.0f; }
static inline void extractPos(const float* m, float& x, float& y, float& z) { x = m[12]; y = m[13]; z = m[14]; }

static std::string resolvePath(const std::string& raw) {
    if (raw.rfind("mdl::", 0) == 0) {
        std::string p = MODELS_ROOT + raw.substr(5);
        std::cout << "[RESOLVE] mdl:: -> " << p << "\n";
        return p;
    }
    if (raw.rfind("msh::", 0) == 0) {
        std::string p = MESHES_ROOT + raw.substr(5);
        std::cout << "[RESOLVE] msh:: -> " << p << "\n";
        return p;
    }
    if (raw.rfind("meshes:", 0) == 0) {
        std::string p = MESHES_ROOT + raw.substr(7);
        std::cout << "[RESOLVE] meshes: -> " << p << "\n";
        return p;
    }
    std::cout << "[RESOLVE] raw -> " << raw << "\n";
    return raw;
}

/*
static bool AppendFromN3WithTransform(const std::string& n3Path,
    const glm::vec3& jsonPos,
    const glm::vec3& jsonRot,
    const glm::vec3& jsonScale)
{
    size_t before = sceneDraws.size();
    std::cout << "[N3] open " << n3Path << "\n";

    Parser parser;
    if (!parser.parse_file(n3Path)) {
        errout << "[N3] parse fail " << n3Path << "\n";
        return false;
    }

    for (auto& n : parser.getNodes()) n->buildLocalMatrix();
    computeHierarchyMatrices(parser.getNodes());

    int added = 0;
    for (auto& node : parser.getNodes()) {
        if (node->mesh_resource_id.empty()) continue;

        DrawCmd dc;
        dc.group = node->primitive_group_idx;

        dc.tex[0] = loadNebTex(node->shader_textures, "DiffMap0");
        dc.tex[1] = loadNebTex(node->shader_textures, "SpecMap0");
        dc.tex[2] = loadNebTex(node->shader_textures, "BumpMap0");
        dc.tex[3] = loadNebTex(node->shader_textures, "EmsvMap0");
        for (int i = 0; i < 4; ++i) {
            dc.has[i] = (dc.tex[i] != 0);
            if (!dc.has[i]) {
                if (i == 0) dc.tex[i] = gWhiteTex;
                else dc.tex[i] = gBlackTex;
            }
        }

        std::string meshPath = node->mesh_resource_id;
        if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
        meshPath = MESHES_ROOT + meshPath;
        if (!LoadNVX2(meshPath, dc.mesh)) {
            errout << "[N3] nvx2 fail " << meshPath << "\n";
            continue;
        }

        for (int s = 0; s < 4; ++s) {
            dc.uvXform[s][0] = 1;
            dc.uvXform[s][1] = 1;
            dc.uvXform[s][2] = 0;
            dc.uvXform[s][3] = 0;
            dc.uvSet[s] = 0;
        }

        glm::vec3 n3Scale(node->scale.x, node->scale.y, node->scale.z);
        glm::vec3 combinedScale(
            jsonScale.x * n3Scale.x,
            jsonScale.y * n3Scale.y,
            jsonScale.z * n3Scale.z
        );
        glm::mat4 finalM = composeTRS(
            jsonPos,
            glm::vec3(node->rotation.x, node->rotation.y, node->rotation.z),
            glm::vec3(node->scale.x, node->scale.y, node->scale.z)
        );

        matCopy(dc.worldMatrix, finalM);

        float px, py, pz;
        extractPos(dc.worldMatrix, px, py, pz);
        errout << "[ADD N3] " << meshPath
            << " pos(" << px << ", " << py << ", " << pz << ")"
            << " jsonScale(" << jsonScale.x << "," << jsonScale.y << "," << jsonScale.z << ")"
            << " n3Scale(" << n3Scale.x << "," << n3Scale.y << "," << n3Scale.z << ")"
            << " final(" << combinedScale.x << "," << combinedScale.y << "," << combinedScale.z << ")\n";

        sceneDraws.push_back(std::move(dc));
        added++;
    }

    std::cout << "[N3] added " << added << " from " << n3Path << "\n";
    return sceneDraws.size() > before;
}
*/

static bool AppendFromNVX2WithTransform(const std::string& nvx2Path, const glm::mat4& rootM) {
    std::cout << "[NVX2] open " << nvx2Path << "\n";

    DrawCmd dc;
    dc.group = -1;

    if (!LoadNVX2(nvx2Path, dc.mesh)) {
        errout << "[NVX2] load fail " << nvx2Path << "\n";
        return false;
    }

    for (int i = 0; i < 4; ++i) {
        dc.uvXform[i][0] = 1;
        dc.uvXform[i][1] = 1;
        dc.uvXform[i][2] = 0;
        dc.uvXform[i][3] = 0;
        dc.uvSet[i] = 0;
    }

    dc.tex[0] = gWhiteTex;
    dc.tex[1] = gBlackTex;
    dc.tex[2] = gFlatNormalTex;
    dc.tex[3] = gBlackTex;
    dc.has[0] = 1; dc.has[1] = 0; dc.has[2] = 0; dc.has[3] = 0;

    matCopy(dc.worldMatrix, rootM);

    float px, py, pz;
    extractPos(dc.worldMatrix, px, py, pz);
    std::cout << "[ADD MESH] " << nvx2Path
        << " pos(" << px << ", " << py << ", " << pz
        << ") Diff=" << dc.tex[0] << " Norm=" << dc.tex[2] << "\n";

    withTransform.push_back(std::move(dc));
    return true;
}


bool LoadSceneFromJson(const std::string& jsonFile) {
    std::ifstream in(jsonFile);
    if (!in) { errout << "[SCENE] open fail " << jsonFile << "\n"; return false; }

    nlohmann::json j;
    try { in >> j; }
    catch (...) { errout << "[SCENE] json parse fail " << jsonFile << "\n"; return false; }
    if (!j.contains("objects") || !j["objects"].is_array()) { errout << "[SCENE] no objects\n"; return false; }

    size_t total = j["objects"].size();
    std::cout << "[SCENE] objects=" << total << "\n";

    int added = 0, skipped = 0, idx = 0;
    for (auto& o : j["objects"]) {
        idx++;
        if (!o.contains("path")) { skipped++; continue; }

        glm::vec3 jsonPos(0.0f), jsonRot(0.0f);
        if (o.contains("pos") && o["pos"].is_array() && o["pos"].size() >= 3) {
            jsonPos = glm::vec3(o["pos"][0].get<float>(), o["pos"][1].get<float>(), o["pos"][2].get<float>());
        }
        if (o.contains("rot") && o["rot"].is_array() && o["rot"].size() >= 3) {
            jsonRot = glm::vec3(o["rot"][0].get<float>(), o["rot"][1].get<float>(), o["rot"][2].get<float>());
        }
        glm::mat4 Mroot = composeTRS(jsonPos, jsonRot, glm::vec3(1.0f));

        std::string raw = o["path"].get<std::string>();
        std::string resolved = resolvePath(raw);

        if (raw.rfind("mdl::", 0) == 0) {
            Parser parser;
            if (!parser.parse_file(resolved)) { skipped++; continue; }

            for (auto& node : parser.getNodes()) node->buildLocalMatrix();
            computeHierarchyMatrices(parser.getNodes());

            for (auto& node : parser.getNodes()) {
                if (node->mesh_resource_id.empty()) continue;

                DrawCmd dc;
                dc.group = node->primitive_group_idx;

                dc.tex[0] = loadNebulaByKeys(node->shader_textures, { "DiffMap0" });
                dc.tex[1] = loadNebulaByKeys(node->shader_textures, { "SpecMap0"});
                dc.tex[2] = loadNebulaByKeys(node->shader_textures, { "BumpMap0" });
                dc.tex[3] = loadNebulaByKeys(node->shader_textures, { "EmsvMap0" });

                for (int i = 0; i < 4; ++i) {
                    dc.has[i] = (dc.tex[i] != 0);
                    if (!dc.has[i]) {
                        if (i == 0) dc.tex[i] = gWhiteTex; 
                        else if (i == 2) dc.tex[i] = gFlatNormalTex;
                        else dc.tex[i] = gBlackTex;
                    }
                }

                std::string meshPath = node->mesh_resource_id;
                if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
                meshPath = MESHES_ROOT + meshPath;
                if (!LoadNVX2(meshPath, dc.mesh)) { errout << "[N3] nvx2 fail " << meshPath << "\n"; continue; }

                for (int s = 0; s < 4; ++s) {
                    dc.uvXform[s][0] = 1; dc.uvXform[s][1] = 1; dc.uvXform[s][2] = 0; dc.uvXform[s][3] = 0;
                    dc.uvSet[s] = 0;
                    auto it = node->uvXformBySlot.find(s);
                    if (it != node->uvXformBySlot.end()) {
                        dc.uvXform[s][0] = it->second[0];
                        dc.uvXform[s][1] = it->second[1];
                        dc.uvXform[s][2] = it->second[2];
                        dc.uvXform[s][3] = it->second[3];
                    }
                    auto iu = node->uvSetBySlot.find(s);
                    if (iu != node->uvSetBySlot.end()) dc.uvSet[s] = iu->second;
                }

                glm::mat4 nodeM = matFromPtr(node->worldMatrix);
                glm::mat4 finalM = Mroot * nodeM;
                matCopy(dc.worldMatrix, finalM);

                withTransform.push_back(std::move(dc));
                added++;
            }
        }
        else if (raw.rfind("msh::", 0) == 0 || raw.rfind("meshes:", 0) == 0) {
            glm::mat4 M = composeTRS(jsonPos, jsonRot, glm::vec3(1.0f));
            if (AppendFromNVX2WithTransform(resolved, M)) added++; else skipped++;
        }
        else {
            skipped++;
        }
    }

    std::cout << "[SCENE] added=" << added << " skipped=" << skipped << " total=" << total << "\n";
    return added > 0;
}


static void DumpSceneDraws() {
   // std::cout << "[DUMP] sceneDraws=" << sceneDraws.size() << "\n";
    for (size_t i = 0; i < withTransform.size(); ++i) {
        const auto& d = withTransform[i];
        float x, y, z; extractPos(d.worldMatrix, x, y, z);
      //  std::cout << "[DRAW] i=" << i << " pos(" << x << ", " << y << ", " << z << ") group=" << d.group << " idxCount=" << (int)d.mesh.idx.size() << " vtxCount=" << (int)d.mesh.verts.size() << "\n";
    }
}

const char* gridVertexShaderSource = R"(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
}
)";

const char* gridFragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
uniform vec3 objectColor;
void main() {
    FragColor = vec4(objectColor, 1.0);
}
)";


class FPSOverlay {
public:
    void Init(int winW, int winH) {
        winWidth = winW;
        winHeight = winH;
        lastTime = glfwGetTime();

        const char* vsSrc = R"(
            #version 330 core
            layout(location=0) in vec2 aPos;
            uniform vec2 uResolution;
            void main() {
                vec2 pos = aPos / uResolution * 2.0 - 1.0;
                gl_Position = vec4(pos.x, -pos.y, 0.0, 1.0);
            }
        )";

        const char* fsSrc = R"(
            #version 330 core
            out vec4 FragColor;
            uniform vec3 uColor;
            void main() {
                FragColor = vec4(uColor, 1.0);
            }
        )";

        shader = CompileShader(vsSrc, fsSrc);

        glGenVertexArrays(1, &vao);
        glGenBuffers(1, &vbo);
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferData(GL_ARRAY_BUFFER, 4096 * sizeof(float), nullptr, GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    void Update() {
        frames++;
        double currentTime = glfwGetTime();
        if (currentTime - lastTime >= 1.0) {
            snprintf(fpsText, sizeof(fpsText), "FPS:%d", frames);
            frames = 0;
            lastTime += 1.0;
        }
    }

    void Render() {
        std::vector<float> verts;
        BuildTextVerts(fpsText, verts);

        glUseProgram(shader);
        glUniform3f(glGetUniformLocation(shader, "uColor"), 1.0f, 1.0f, 0.0f);
        glUniform2f(glGetUniformLocation(shader, "uResolution"), (float)winWidth, (float)winHeight);

        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glBufferSubData(GL_ARRAY_BUFFER, 0, verts.size() * sizeof(float), verts.data());
        glDrawArrays(GL_TRIANGLES, 0, verts.size() / 2);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    void Cleanup() {
        glDeleteProgram(shader);
        glDeleteVertexArrays(1, &vao);
        glDeleteBuffers(1, &vbo);
    }

private:
    GLuint vao = 0, vbo = 0, shader = 0;
    int winWidth = 0, winHeight = 0;
    int frames = 0;
    double lastTime = 0.0;
    char fpsText[32] = "FPS:0";

    GLuint CompileShader(const char* vsSrc, const char* fsSrc) {
        GLuint vs = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vs, 1, &vsSrc, nullptr);
        glCompileShader(vs);

        GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fs, 1, &fsSrc, nullptr);
        glCompileShader(fs);

        GLuint prog = glCreateProgram();
        glAttachShader(prog, vs);
        glAttachShader(prog, fs);
        glLinkProgram(prog);

        glDeleteShader(vs);
        glDeleteShader(fs);
        return prog;
    }

    // glyph bitmaps
    static const unsigned char font_digits[10][8];
    static const unsigned char font_F[8];
    static const unsigned char font_P[8];
    static const unsigned char font_S[8];
    static const unsigned char font_colon[8];

    void BuildGlyphVerts(const unsigned char glyph[8], float x, float y, float scale, std::vector<float>& verts) {
        for (int row = 0; row < 8; ++row) {
            unsigned char bits = glyph[row];
            for (int col = 0; col < 8; ++col) {
                if (bits & (1 << (7 - col))) {
                    float px = x + col * scale;
                    float py = y + row * scale;
                    float s = scale;
                    verts.insert(verts.end(), { px,py,  px + s,py,  px + s,py + s });
                    verts.insert(verts.end(), { px,py,  px + s,py + s,  px,py + s });
                }
            }
        }
    }

    void BuildTextVerts(const char* text, std::vector<float>& verts) {
        float x = 10, y = 10, scale = 2.0f;
        for (const char* p = text; *p; ++p) {
            char c = *p;
            if (c >= '0' && c <= '9') {
                BuildGlyphVerts(font_digits[c - '0'], x, y, scale, verts);
                x += 8 * scale;
            }
            else if (c == 'F') { BuildGlyphVerts(font_F, x, y, scale, verts); x += 8 * scale; }
            else if (c == 'P') { BuildGlyphVerts(font_P, x, y, scale, verts); x += 8 * scale; }
            else if (c == 'S') { BuildGlyphVerts(font_S, x, y, scale, verts); x += 8 * scale; }
            else if (c == ':') {
                x -= 3 * scale;
                BuildGlyphVerts(font_colon, x, y, scale, verts);
                x += 8 * scale;
            }
        }
    }
};

// glyph definitions
const unsigned char FPSOverlay::font_digits[10][8] = {
    {0x3C,0x66,0x6E,0x7E,0x76,0x66,0x3C,0x00},
    {0x18,0x38,0x18,0x18,0x18,0x18,0x7E,0x00},
    {0x3C,0x66,0x06,0x1C,0x30,0x66,0x7E,0x00},
    {0x3C,0x66,0x06,0x1C,0x06,0x66,0x3C,0x00},
    {0x0C,0x1C,0x3C,0x6C,0x7E,0x0C,0x0C,0x00},
    {0x7E,0x60,0x7C,0x06,0x06,0x66,0x3C,0x00},
    {0x1C,0x30,0x60,0x7C,0x66,0x66,0x3C,0x00},
    {0x7E,0x66,0x0C,0x18,0x18,0x18,0x18,0x00},
    {0x3C,0x66,0x66,0x3C,0x66,0x66,0x3C,0x00},
    {0x3C,0x66,0x66,0x3E,0x06,0x0C,0x38,0x00}
};
const unsigned char FPSOverlay::font_F[8] = { 0x7E,0x60,0x60,0x7C,0x60,0x60,0x60,0x00 };
const unsigned char FPSOverlay::font_P[8] = { 0x7C,0x66,0x66,0x7C,0x60,0x60,0x60,0x00 };
const unsigned char FPSOverlay::font_S[8] = { 0x3C,0x60,0x60,0x3C,0x06,0x06,0x7C,0x00 };
const unsigned char FPSOverlay::font_colon[8] = { 0x00,0x18,0x18,0x00,0x18,0x18,0x00,0x00 };

std::vector<DrawCmd> BuildDrawsWithTransform(const Parser& parser,
    const glm::vec3& pos,
    const glm::vec3& rot,
    const glm::vec3& scale) {
    std::vector<DrawCmd> out;
    glm::mat4 trs = composeTRS(pos, rot, scale);

    for (auto& node : parser.getNodes()) {
        if (node->mesh_resource_id.empty()) {
            continue;
        }

        DrawCmd dc;
        dc.group = node->primitive_group_idx;

        glm::mat4 nodeLocal;
        buildTransformMatrix(&nodeLocal[0][0],
            { node->position.x, node->position.y, node->position.z },
            { node->rotation.x, node->rotation.y, node->rotation.z, node->rotation.w },
            { node->scale.x, node->scale.y, node->scale.z });
        glm::mat4 finalM = trs * nodeLocal;
        matCopy(dc.worldMatrix, finalM);

        dc.tex[0] = loadNebTex(node->shader_textures, "DiffMap0");
        dc.tex[1] = loadNebTex(node->shader_textures, "SpecMap0");
        dc.tex[2] = loadNebTex(node->shader_textures, "BumpMap0");
        dc.tex[3] = loadNebTex(node->shader_textures, "EmsvMap0");
        for (int i = 0; i < 4; ++i) {
            dc.has[i] = (dc.tex[i] != 0);
            if (!dc.has[i]) dc.tex[i] = (i == 0) ? gWhiteTex : gBlackTex;
        }

        std::string meshPath = node->mesh_resource_id;
        if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
        meshPath = MESHES_ROOT + meshPath;
        if (!LoadNVX2(meshPath, dc.mesh)) {
            continue;
        }

        dc.name = meshPath;

        for (int s = 0; s < 4; ++s) {
            dc.uvXform[s][0] = 1;
            dc.uvXform[s][1] = 1;
            dc.uvXform[s][2] = 0;
            dc.uvXform[s][3] = 0;
            dc.uvSet[s] = 0;
        }
        out.push_back(std::move(dc));
    }
    return out;
}

static inline std::string strlower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return s;
}
 
static bool LoadNVX2IntoScene(const std::string& path) {
    glm::mat4 I(1.0f);
    ensureFallbacks();
    if (!AppendFromNVX2WithTransform(path, I)) {
        std::cout << "[DROP] NVX2 load failed: " << path << "\n";
        return false;
    }
    if (!withTransform.empty()) {
        auto& dc = withTransform.back();
        for (int i = 0; i < 4; i++) {
            dc.tex[i] = gWhiteTex;
            dc.has[i] = 1;
        }
        std::cout << "[DROP] NVX2 forced white textures: " << path
            << " texIDs(" << dc.tex[0] << "," << dc.tex[1]
            << "," << dc.tex[2] << "," << dc.tex[3] << ")\n";
    }
    else {
        std::cout << "[DROP] NVX2 loaded but withTransform empty: " << path << "\n";
    }
    return true;
}


void appendN3WTransform(const std::string& path,
    const glm::vec3& pos,
    const glm::vec3& rot,
    const glm::vec3& scale)
{
    Parser parser;
    if (!parser.parse_file(path)) return;
    for (auto& node : parser.getNodes()) node->buildLocalMatrix();
    computeHierarchyMatrices(parser.getNodes());

    auto draws = BuildDrawsWithTransform(parser, pos, rot, scale);
    for (auto& dc : draws) withTransform.push_back(std::move(dc));
}

void InitDDSQuad() {
    float size = 0.5f;
    float quadVertices[] = {
        -size,  size,  0.0f, 1.0f,
        -size, -size,  0.0f, 0.0f,
         size, -size,  1.0f, 0.0f,

        -size,  size,  0.0f, 1.0f,
         size, -size,  1.0f, 0.0f,
         size,  size,  1.0f, 1.0f
    };
    glGenVertexArrays(1, &gDDSQuadVAO);
    glGenBuffers(1, &gDDSQuadVBO);
    glBindVertexArray(gDDSQuadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, gDDSQuadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);

}

void RenderDDSQuad() {
    if (!gShowDDS || gDDSTexture == 0) return;
    glDisable(GL_DEPTH_TEST);
    glUseProgram(gDDSShader);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, gDDSTexture);
    glUniform1i(glGetUniformLocation(gDDSShader, "tex"), 0);
    glBindVertexArray(gDDSQuadVAO);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);
    glEnable(GL_DEPTH_TEST);
}

static void HandleDroppedFile(const std::string& path, bool& didClear) {
    std::string l = strlower(path);
    if (l.size() >= 5 && l.compare(l.size() - 5, 5, ".nvx2") == 0) {
        if (!didClear) { withTransform.clear(); didClear = true; }
        LoadNVX2IntoScene(path);
        return;
    }
    if (l.size() >= 3 && l.compare(l.size() - 3, 3, ".n3") == 0) {
        if (!didClear) { withTransform.clear(); didClear = true; }
        appendN3WTransform(path, { 0,0,0 }, { 0,0,0 }, { 1,1,1 });
        return;
    }
    if (l.size() >= 4 && l.compare(l.size() - 4, 4, ".dds") == 0) {
        unsigned int tex = LoadDDS(path);
        if (tex != 0) {
            gDDSTexture = tex;
            if (!gShowDDS) {
                logout << "Quad inactive, no textures visible!\n";
            }
            else {
                std::cout << "[DROP] DDS fullscreen quad: " << path << " id=" << tex << "\n";
            }
        }
        return;
    }
    std::cout << "[DROP] Unsupported: " << path << "\n";
}

static void DropCallback(GLFWwindow* window, int count, const char** paths) {
    bool didClear = false;
    for (int i = 0; i < count; i++) HandleDroppedFile(paths[i], didClear);
}

static void InstallDragDrop(GLFWwindow* window) {
    glfwSetDropCallback(window, DropCallback);
}


static GLuint makeSolidTex(unsigned char r, unsigned char g, unsigned char b) {
    GLuint t; glGenTextures(1, &t);
    glBindTexture(GL_TEXTURE_2D, t);
    unsigned char px[3] = { r,g,b };
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB8, 1, 1, 0, GL_RGB, GL_UNSIGNED_BYTE, px);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    return t;
}

void ensureFallbacks() {
    if (!gWhiteTex) {
        uint32_t whitePixel = 0xFFFFFFFF;
        glGenTextures(1, &gWhiteTex);
        glBindTexture(GL_TEXTURE_2D, gWhiteTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0,
            GL_RGBA, GL_UNSIGNED_BYTE, &whitePixel);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    if (!gBlackTex) {
        uint32_t blackPixel = 0xFF000000;
        glGenTextures(1, &gBlackTex);
        glBindTexture(GL_TEXTURE_2D, gBlackTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0,
            GL_RGBA, GL_UNSIGNED_BYTE, &blackPixel);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
    if (!gFlatNormalTex) {
        uint32_t flatNormalPixel = 0xFF8080FF;
        glGenTextures(1, &gFlatNormalTex);
        glBindTexture(GL_TEXTURE_2D, gFlatNormalTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0,
            GL_RGBA, GL_UNSIGNED_BYTE, &flatNormalPixel);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glBindTexture(GL_TEXTURE_2D, 0);
    }
}

GLuint loadNebulaByKeys(const std::unordered_map<std::string, std::string>& m,
    std::initializer_list<const char*> keys) {
    for (auto k : keys) {
        auto it = m.find(k);
        if (it == m.end()) continue;
        std::string full = resolveTexture(it->second);
        auto itc = gTexCache.find(full);
        if (itc != gTexCache.end()) return itc->second;
        GLuint t = LoadDDS(full);
        if (t) {
            gTexCache[full] = t;
            if (debug_output) std::cout << "[TEX] " << k << " -> " << t << " (" << full << ")\n";
            return t;
        }
    }
    return 0;
}

struct InstanceTRS { glm::vec3 pos, rotDeg, scale; };



struct OverrideTRS {
    glm::vec3 pos;
    glm::vec3 rot;
    glm::vec3 scale;
};
static std::unordered_map<std::string, OverrideTRS> gOverrides = {
    
    { "deco_rubble_bruned_01.n3",     { { 2.220055f, 1.986336f, 2.411935f }, { 104.07255f, 80.529389f, 79.495329f }, { 1.039876442f, 1.00074181176f, 0.99903148f } } },
    { "fx_fire_02_small.n3",          { { 1.000000f, 1.000000f, 1.000000f }, { -102.25024f, 14.796412f, 78.501279f }, { 0.0028649580f, 0.077471413f, 0.011743932f } } },

    { "tile_catwalk_fence_01.n3",     { { 0.379757f, 4.234385f, 2.021405f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "tile_catwalk_fence_02.n3",     { { 0.357390f, 4.265504f, 1.971051f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },

    { "char_screen_stone_streaked_01.n3", { { 0.655783f, 0.314703f, 0.510550f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "char_screen_stone_02.n3",      { { 0.243499f, 0.226329f, 0.228052f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "char_screen_stone_01.n3",      { { 0.199286f, 0.206600f, 0.213412f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },

    { "deco_stone_plate_01.n3",       { { 1.500604f, 0.527159f, 2.023414f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "deco_stone_03.n3",             { { 0.662257f, 0.635154f, 0.933427f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },

    { "veg_root_01.n3",               { { 2.356448f, 0.844028f, 2.493979f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "veg_grass_03.n3",              { { 2.321698f, 0.611534f, 2.474809f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },

    { "fx_neferati_fog_04.n3",        { { 3.504755f, 1.000000f, 2.187984f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },
    { "fx_nebulapatch_01.n3",         { { 3.947490f, 0.000101f, 3.258760f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } },

    { "stone_platish_01.n3",          { { 0.551777f, 0.197556f, 0.408828f }, { 0.0f, 0.0f, 0.0f }, { 1.0f, 1.0f, 1.0f } } }
};



static void ApplyOverride(const std::string& path, DrawCmd& dc) {
    for (auto& [key, trs] : gOverrides) {
        if (path.find(key) != std::string::npos) {
            glm::mat4 m(1.0f);
            m = glm::translate(m, trs.pos);
            m = glm::rotate(m, glm::radians(trs.rot.x), glm::vec3(1, 0, 0));
            m = glm::rotate(m, glm::radians(trs.rot.y), glm::vec3(0, 1, 0));
            m = glm::rotate(m, glm::radians(trs.rot.z), glm::vec3(0, 0, 1));
            m = glm::scale(m, trs.scale);
            matCopy(dc.worldMatrix, m);
            break;
        }
    }
}

std::vector<DrawCmd> BuildDrawsFromNodes(const Parser& parser) {
    std::vector<DrawCmd> out;
    for (auto& node : parser.getNodes()) {
        if (node->mesh_resource_id.empty())
        {
            continue;
        }

        DrawCmd dc;
        dc.group = node->primitive_group_idx;
        memcpy(dc.worldMatrix, node->worldMatrix, sizeof(dc.worldMatrix));

        dc.tex[0] = loadNebTex(node->shader_textures, "DiffMap0");
        dc.tex[1] = loadNebTex(node->shader_textures, "SpecMap0");
        dc.tex[2] = loadNebTex(node->shader_textures, "BumpMap0");
        dc.tex[3] = loadNebTex(node->shader_textures, "EmsvMap0");
        for (int i = 0; i < 4; ++i) {
            dc.has[i] = (dc.tex[i] != 0);
            if (!dc.has[i]) dc.tex[i] = (i == 0) ? gWhiteTex : gBlackTex;
        }

        std::string meshPath = node->mesh_resource_id;
        if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
        meshPath = MESHES_ROOT + meshPath;
        if (!LoadNVX2(meshPath, dc.mesh)) {

            continue;
        }

        dc.name = meshPath;

        for (int s = 0; s < 4; ++s) {
            auto it = node->uvXformBySlot.find(s);
            if (it != node->uvXformBySlot.end()) {
                dc.uvXform[s][0] = it->second[0];
                dc.uvXform[s][1] = it->second[1];
                dc.uvXform[s][2] = it->second[2];
                dc.uvXform[s][3] = it->second[3];
            }
            else {
                dc.uvXform[s][0] = 1;
                dc.uvXform[s][1] = 1;
                dc.uvXform[s][2] = 0;
                dc.uvXform[s][3] = 0;
            }
            auto iu = node->uvSetBySlot.find(s);
            dc.uvSet[s] = (iu != node->uvSetBySlot.end()) ? iu->second : 0;
        }
        out.push_back(std::move(dc));
    }
    return out;
}



static const unsigned char font8x8_basic[128][8] = {
    /* 48 '0' */
    {0x3C,0x66,0x6E,0x7E,0x76,0x66,0x3C,0x00},
    /* 49 '1' */
    {0x18,0x38,0x18,0x18,0x18,0x18,0x7E,0x00},
    /* 50 '2' */
    {0x3C,0x66,0x06,0x1C,0x30,0x66,0x7E,0x00},
    /* 51 '3' */
    {0x3C,0x66,0x06,0x1C,0x06,0x66,0x3C,0x00},
    /* 52 '4' */
    {0x0C,0x1C,0x3C,0x6C,0x7E,0x0C,0x0C,0x00},
    /* 53 '5' */
    {0x7E,0x60,0x7C,0x06,0x06,0x66,0x3C,0x00},
    /* 54 '6' */
    {0x1C,0x30,0x60,0x7C,0x66,0x66,0x3C,0x00},
    /* 55 '7' */
    {0x7E,0x66,0x0C,0x18,0x18,0x18,0x18,0x00},
    /* 56 '8' */
    {0x3C,0x66,0x66,0x3C,0x66,0x66,0x3C,0x00},
    /* 57 '9' */
    {0x3C,0x66,0x66,0x3E,0x06,0x0C,0x38,0x00},
    /* 58 ':' */
    {0x00,0x18,0x18,0x00,0x18,0x18,0x00,0x00},
    /* 70 'F' */
    {0x7E,0x60,0x60,0x7C,0x60,0x60,0x60,0x00},
    /* 80 'P' */
    {0x7C,0x66,0x66,0x7C,0x60,0x60,0x60,0x00},
    /* 83 'S' */
    {0x3C,0x60,0x60,0x3C,0x06,0x06,0x7C,0x00},
};


void drawChar8x8(float x, float y, char c) {
    if (c < 0 || c > 127) return;
    glBegin(GL_QUADS);
    for (int row = 0; row < 8; row++) {
        unsigned char bits = font8x8_basic[(int)c][row];
        for (int col = 0; col < 8; col++) {
            if (bits & (1 << (7 - col))) {
                float px = x + col * 2;
                float py = y + row * 2;
                glVertex2f(px, py);
                glVertex2f(px + 2, py);
                glVertex2f(px + 2, py + 2);
                glVertex2f(px, py + 2);
            }
        }
    }
    glEnd();
}

void drawText8x8(float x, float y, const char* text) {
    for (const char* p = text; *p; p++, x += 10) {
        drawChar8x8(x, y, *p);
    }
}

const char* quadDDSVert = R"(
#version 330 core
layout (location = 0) in vec2 aPos;
layout (location = 1) in vec2 aUV;
out vec2 vUV;
void main(){
    vUV = aUV;
    gl_Position = vec4(aPos, 0.0, 1.0);
}
)";

const char* quadDDSFrag = R"(
#version 330 core
in vec2 vUV;
out vec4 FragColor;
uniform sampler2D tex;
void main(){
    FragColor = texture(tex, vUV);
}
)";


int main(int argc, char* argv[])
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_ANY_PROFILE);

    double sizeX = 1200, sizeY = 800;
    GLFWwindow* window = glfwCreateWindow(static_cast<int>(sizeX), static_cast<int>(sizeY), "SIMO'S NEB3 VIEWER", nullptr, nullptr);
    if (!window) { glfwTerminate(); return -1; }

    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, [](GLFWwindow* w, int key, int sc, int action, int mods) {
        if (key >= 0 && key < 1024) { if (action == GLFW_PRESS) keys[key] = true; else if (action == GLFW_RELEASE) keys[key] = false; }
        });
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) return -1;
    glViewport(0, 0, (int)sizeX, (int)sizeY);

    InstallDragDrop(window);

    unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    CheckShader(vertexShader, "VS");

    unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    CheckShader(fragmentShader, "FS");

    unsigned int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    CheckProgram(shaderProgram);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    GLuint gridVS = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(gridVS, 1, &gridVertexShaderSource, NULL);
    glCompileShader(gridVS);
    CheckShader(gridVS, "GridVS");

    GLuint gridFS = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(gridFS, 1, &gridFragmentShaderSource, NULL);
    glCompileShader(gridFS);
    CheckShader(gridFS, "GridFS");

    GLuint gridShaderProgram = glCreateProgram();
    glAttachShader(gridShaderProgram, gridVS);
    glAttachShader(gridShaderProgram, gridFS);
    glLinkProgram(gridShaderProgram);
    CheckProgram(gridShaderProgram);
    glDeleteShader(gridVS);
    glDeleteShader(gridFS);

    GLuint ddsVS = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(ddsVS, 1, &quadDDSVert, NULL);
    glCompileShader(ddsVS);
    CheckShader(ddsVS, "DDSVS");

    GLuint ddsFS = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(ddsFS, 1, &quadDDSFrag, NULL);
    glCompileShader(ddsFS);
    CheckShader(ddsFS, "DDSFS");

    GLuint ddsShaderProgram = glCreateProgram();
    glAttachShader(ddsShaderProgram, ddsVS);
    glAttachShader(ddsShaderProgram, ddsFS);
    glLinkProgram(ddsShaderProgram);
    CheckProgram(ddsShaderProgram);
    glDeleteShader(ddsVS);
    glDeleteShader(ddsFS);

    gDDSShader = ddsShaderProgram;


    InitGrid(10.0f, 1.0f);

    std::vector<std::string> n3Files = {
        /*
        "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_04.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_05.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_rock_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_rock_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_rock_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_rock_flat_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_rock_flat_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_stones_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_stones_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_deco_stones_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_destroy_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_01_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_02_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_02_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_03_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_04.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_04_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_fall_corner_04_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_01_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_01_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_02_end.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_liquid_lava_straight_02_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_river_fall_corner_01_ne.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_river_fall_corner_02_ne.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_river_fall_straight_02_ne.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_river_fall_straight_02_nw.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_river_fall_straight_02_sw.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_floor_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_floor_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_bridge_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_curve_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_curve_01_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_curve_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_curve_02_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_end_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_end_01_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_01_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_01_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_02_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_02_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_03_ext_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_fall_03_lava.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_ramp_corner_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_ramp_corner_05.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_ramp_straight_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_wall_connection.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_wall_connection_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_wall_corner_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_01_var_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_02_var_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_04.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_05.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_bush_06.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_04.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_05.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_tree_pine_06.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_01.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_02.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_mushroom_03.n3",
        "C:/drasa_online/work/models/m006_imperial_lava/m006_veg_zypresse_01.n3",


       // "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_03.n3",
      // "C:/drasa_online/work/models/m006_imperial_lava/m006_decal_ground_05.n3",
       // "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_wall_straight_03_var_02.n3",
       // "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_lava_straight_03_lava.n3",
      //  "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ramp_rock_02_ext_01.n3",
      //  "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_floor_03.n3",
      //  "C:/drasa_online/work/models/m006_imperial_lava/m006_tile_ground_wall_rock_03.n3",

        */

        "C:/drasa_online/work/models/t000_char/background_01.n3",  
        "C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        "C:/drasa_online/work/models/t000_char/char_screen_stone_02.n3",
        "C:/drasa_online/work/models/t000_char/char_screen_stone_streaked_01.n3",
        "C:/drasa_online/work/models/t000_char/decal_sand_01.n3",
        "C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        "C:/drasa_online/work/models/t000_char/deco_stone_plate_01.n3",
         "C:/drasa_online/work/models/m016_mountain_fortress/deco_rubble_bruned_01.n3",
        "C:/drasa_online/work/models/t000_char/fx_lightbeam_01.n3",
        "C:/drasa_online/work/models/t000_char/fx_nebulapatch_01.n3",
        "C:/drasa_online/work/models/t000_char/ground_02.n3",
        "C:/drasa_online/work/models/t000_char/stone_platish_01.n3",
        "C:/drasa_online/work/models/t000_char/veg_grass_03.n3",
        "C:/drasa_online/work/models/t000_char/veg_root_01.n3",
        "C:/drasa_online/work/models/t001_hub/fx_lightbeam_01.n3",
        "C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        "C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
       "C:/drasa_online/work/models/t010_caribbean_01/fx_neferati_fog_04.n3"
    };

    /*
    for (auto& n3FullPath : n3Files) {
        Parser parser;
        if (!parser.parse_file(n3FullPath)){
            errout << "Cannot parse\n";
            continue;
        }
        for (auto& node : parser.getNodes()) node->buildLocalMatrix();
        computeHierarchyMatrices(parser.getNodes());

        auto draws = BuildDrawsFromNodes(parser);
        for (auto& dc : draws) {
           // ApplyOverride(n3FullPath, dc);
           // fromNodes.push_back(std::move(dc));
        }
    }*/

    auto appendN3WTransform = [&](const std::string& path,
        const glm::vec3& pos,
        const glm::vec3& rot,
        const glm::vec3& scale)
        {
            Parser parser;
            if (!parser.parse_file(path)) return;
            for (auto& node : parser.getNodes()) node->buildLocalMatrix();
            computeHierarchyMatrices(parser.getNodes());

            auto draws = BuildDrawsWithTransform(parser, pos, rot, scale);
            for (auto& dc : draws) withTransform.push_back(std::move(dc));
        };

    //appendN3WTransform("C:/drasa_online/work/models/t000_char/ground_02.n3",
    //    { 0.207103,0.125945,-1.487267 }, { 0,0,0 }, { 1,1,1 });

    /*
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/arch_house_01.n3",
        { -5,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/arch_house_01.n3",
        { 0,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/arch_house_01.n3",
        { 5,0,0 }, { 0,0,0 }, { 1,1,1 });
    
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_ground_stone_01_var_01.n3",
        { 0,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_ground_stone_border_01.n3",
        { 5,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_ground_stone_corner_01.n3",
        { 10,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_pier_bar_01.n3",
        { 15,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_pier_bar_02.n3",
        { 20,0,0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_ramp_03.n3",
        { 25,0,0 }, { 0,0,0 }, { 1,1,1 });

    appendN3WTransform("C:/drasa_online/work/models/m016_mountain_fortress/deco_rubble_bruned_01.n3",
        { 0.207103,0.125945,-1.487267 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -4.929287,-0.062122,-10.496880 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 7.742453,0.606517,-14.063550 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -6.172833,0.694790,-2.424380 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 2.273689,-0.121829,-13.564480 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -5.296444,0.583664,-5.641184 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 3.815803,0.492663,-14.623500 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -6.011575,1.295615,-12.791900 }, { 0,0,0 }, { 1,1,1 });

    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_streaked_01.n3",
        { 3.815803,0.492663,-14.623500 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { -6.011575,1.295615,-12.791900 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { 6.787735,-0.046394,-10.666000 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_plate_01.n3",
        { 5.115231,-0.174960,-10.977040 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/stone_platish_01.n3",
        { 5.618426,0.946000,-11.598300 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_streaked_01.n3",
        { 5.565071,0.159313,-4.693389 }, { 0,0,0 }, { 1,1,1 });

    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 1.319173,-0.076579,-7.056397 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 1.319173,-0.076579,-7.056397 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -1.839560,-0.049209,-0.832065 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -0.504863,-0.585045,-15.265860 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -1.333145,-0.670562,1.624351 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 2.357347,-0.049209,-1.453099 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 2.357347,-0.049209,-4.259654 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { 1.061164,-1.507582,3.526766 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/char_screen_stone_01.n3",
        { -2.522160,-0.049209,-4.257953 }, { 0,0,0 }, { 1,1,1 });

    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { -4.502991,1.460832,-5.866385 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { 5.320754,1.779095,-4.419336 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { 3.009094,0.651966,-12.839520 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { 2.665146,-0.066272,-15.897940 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t000_char/deco_stone_03.n3",
        { 3.771819,4.060757,-19.829241 }, { 0,0,0 }, { 1,1,1 });


    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.95, -7.625 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.625, -6.0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.625, -4.375 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.625, -2.75 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.625, -1.125 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -3.625, 0.5 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { 4.75, -0.346, -15.125 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
        { -4.0, -0.346, -7.625 }, { 0,0,0 }, { 1,1,1 });

    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
        { 5.0, 4.318, -7.25 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
        { 5.125, 2.75, -2.875 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
        { 5.625, 3.0, -7.0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
        { -4.563, 2.75, -2.375 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_02.n3",
        { -4.563, 4.0, -19.0 }, { 0,0,0 }, { 1,1,1 });

        */

  //  appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
  //      { -4.976412,-0.039607,0.503584 }, { 0,0,0 }, { 1,1,1 });
  //  appendN3WTransform("C:/drasa_online/work/models/t001_hub/tile_catwalk_fence_01.n3",
  //      { 4.455561,-0.213331,-9.086962 }, { 0,0,0 }, { 1,1,1 });

    logout << "=== fromNodes === Count " << fromNodes.size() << "\n";
    for (auto& dc : fromNodes) {
        float px, py, pz;
        extractPos(dc.worldMatrix, px, py, pz);
       // std::cout << dc.name << " pos(" << px << "," << py << "," << pz << ")\n";
    }

    logout << "\n=== wTransform === Count " << withTransform.size() << "\n";
    for (auto& dc : withTransform) {
        float px, py, pz;
        extractPos(dc.worldMatrix, px, py, pz);
       // std::cout << dc.name << " pos(" << px << "," << py << "," << pz << ")\n";
    }

  //  sceneDraws.clear();
  //  if (!LoadSceneFromJson("C:/drasa_online/work/a0000_char.json")) errout << "[SCENE] load failed\n";
  //  DumpSceneDraws();

    glUseProgram(shaderProgram);

    glUniform1i(glGetUniformLocation(shaderProgram, "DiffMap"), 0);
    glUniform1i(glGetUniformLocation(shaderProgram, "SpecMap"), 1);
    glUniform1i(glGetUniformLocation(shaderProgram, "BumpMap"), 2);
    glUniform1i(glGetUniformLocation(shaderProgram, "EmissiveMap"), 3);

    glUniform3f(glGetUniformLocation(shaderProgram, "SunDir"), -0.3f, -1.0f, -0.2f);
    glUniform3f(glGetUniformLocation(shaderProgram, "SunColor"), 1.0f, 0.95f, 0.85f);
    glUniform1f(glGetUniformLocation(shaderProgram, "SunIntensity"), 2.0f);
    glUniform3f(glGetUniformLocation(shaderProgram, "AmbientColor"), 0.25f, 0.25f, 0.25f);

    glUniform1f(glGetUniformLocation(shaderProgram, "SpecularPower"), 64.0f);
    glUniform1f(glGetUniformLocation(shaderProgram, "SpecularIntensity"), 1.0f);
    glUniform1f(glGetUniformLocation(shaderProgram, "EmissiveIntensity"), 1.0f);

   // glUniform1i(glGetUniformLocation(shaderProgram, "HasBump"), 1);
   // glUniform1i(glGetUniformLocation(shaderProgram, "HasSpec"), 1);

    int uvSet[4] = { 0, 0, 0, 0 };
    glUniform1iv(glGetUniformLocation(shaderProgram, "UvSet"), 4, uvSet);

    float uvXform[16] = {
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0,
        1, 1, 0, 0
    };
    glUniform4fv(glGetUniformLocation(shaderProgram, "UvXform"), 4, uvXform);

    glUseProgram(0);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    std::vector<float> boneId(128 * 16, 0.0f);
    for (int i = 0; i < 128; ++i) {
        boneId[i * 16 + 0] = 1.0f;
        boneId[i * 16 + 5] = 1.0f;
        boneId[i * 16 + 10] = 1.0f;
        boneId[i * 16 + 15] = 1.0f;
    }

    // glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    // glEnable(GL_BLEND);
    // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // glEnable(GL_ALPHA_TEST);

    FPSOverlay fps;
    InitDDSQuad();

    fps.Init((int)sizeX, (int)sizeY);
    while (!glfwWindowShouldClose(window))
    {
        processInput(window);
        fps.Update();

        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);

        float projectionMatrix[16];
        createPerspectiveMatrix(projectionMatrix, 45.0f, (float)sizeX / (float)sizeY, 0.1f, 1500.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, projectionMatrix);

        float viewMatrix[16];
        float camX = sinf(yaw) * cosf(pitch) * distance + targetX;
        float camY = sinf(pitch) * distance + targetY;
        float camZ = cosf(yaw) * cosf(pitch) * distance + targetZ;
        createLookAtMatrix(viewMatrix, camX, camY, camZ, targetX, targetY, targetZ, 0.0f, 1.0f, 0.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, viewMatrix);
        glUniform3f(glGetUniformLocation(shaderProgram, "CameraPos"), camX, camY, camZ);
        glUniform3f(glGetUniformLocation(shaderProgram, "PointPos"), camX, camY, camZ);

        GLint boneLoc = glGetUniformLocation(shaderProgram, "BoneMatrices[0]");
        if (boneLoc >= 0) {
            glUniformMatrix4fv(boneLoc, 128, GL_FALSE, boneId.data());
        }

        // swap withTransform and fromNodes
        for (auto& obj : withTransform) {
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, obj.worldMatrix);

            GLuint d = obj.has[0] ? obj.tex[0] : gWhiteTex;
            GLuint s = obj.has[1] ? obj.tex[1] : gBlackTex;
            GLuint n = obj.has[2] ? obj.tex[2] : gFlatNormalTex;
            GLuint e = obj.has[3] ? obj.tex[3] : gBlackTex;

            glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, d);
            glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, s);
            glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, n);
            glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, e);

            glUniform1i(glGetUniformLocation(shaderProgram, "HasDiff"), obj.has[0] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasSpec"), obj.has[1] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasBump"), obj.has[2] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasEmissive"), obj.has[3] ? 1 : 0);

            glUniform4fv(glGetUniformLocation(shaderProgram, "UvXform[0]"), 4, &obj.uvXform[0][0]);
            glUniform1iv(glGetUniformLocation(shaderProgram, "UvSet[0]"), 4, obj.uvSet);

            obj.applyUVTransforms(shaderProgram);

            glBindVertexArray(obj.mesh.vao);

            if (!obj.mesh.idx.empty()) {
                if (obj.group >= 0 && obj.group < (int)obj.mesh.groups.size()) {
                    const Nvx2Group& g = obj.mesh.groups[obj.group];
                    const GLsizei count = (GLsizei)g.indexCount();
                    const GLsizei offset = (GLsizei)(g.firstIndex() * sizeof(uint32_t));
                    const GLint   base = (GLint)g.baseVertex();
                    glDrawElementsBaseVertex(GL_TRIANGLES, count, GL_UNSIGNED_INT, (void*)(intptr_t)offset, base);
                }
                else {
                    glDrawElements(GL_TRIANGLES, (GLsizei)obj.mesh.idx.size(), GL_UNSIGNED_INT, 0);
                }
            }

            glBindVertexArray(0);
        }

        /*
        for (auto& obj : fromNodes) {
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, obj.worldMatrix);

            GLuint d = obj.has[0] ? obj.tex[0] : gWhiteTex;
            GLuint s = obj.has[1] ? obj.tex[1] : gBlackTex;
            GLuint n = obj.has[2] ? obj.tex[2] : gFlatNormalTex;
            GLuint e = obj.has[3] ? obj.tex[3] : gBlackTex;

            glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, d);
            glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, s);
            glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, n);
            glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, e);

            glUniform1i(glGetUniformLocation(shaderProgram, "HasDiff"), obj.has[0] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasSpec"), obj.has[1] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasBump"), obj.has[2] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasEmissive"), obj.has[3] ? 1 : 0);

            glUniform4fv(glGetUniformLocation(shaderProgram, "UvXform[0]"), 4, &obj.uvXform[0][0]);
            glUniform1iv(glGetUniformLocation(shaderProgram, "UvSet[0]"), 4, obj.uvSet);

            obj.applyUVTransforms(shaderProgram);

            glBindVertexArray(obj.mesh.vao);

            if (!obj.mesh.idx.empty()) {
                if (obj.group >= 0 && obj.group < (int)obj.mesh.groups.size()) {
                    const Nvx2Group& g = obj.mesh.groups[obj.group];
                    const GLsizei count = (GLsizei)g.indexCount();
                    const GLsizei offset = (GLsizei)(g.firstIndex() * sizeof(uint32_t));
                    const GLint   base = (GLint)g.baseVertex();
                    glDrawElementsBaseVertex(GL_TRIANGLES, count, GL_UNSIGNED_INT, (void*)(intptr_t)offset, base);
                }
                else {
                    glDrawElements(GL_TRIANGLES, (GLsizei)obj.mesh.idx.size(), GL_UNSIGNED_INT, 0);
                }
            }

            glBindVertexArray(0);
        }
        */

        DrawGrid(gridShaderProgram, viewMatrix, projectionMatrix);

        fps.Render();

        RenderDDSQuad();

        glUseProgram(0);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    for (auto& obj : withTransform) {
        glDeleteVertexArrays(1, &obj.mesh.vao);
        glDeleteBuffers(1, &obj.mesh.vbo);
        glDeleteBuffers(1, &obj.mesh.ebo);
    }
    for (auto& obj : fromNodes) {
        glDeleteVertexArrays(1, &obj.mesh.vao);
        glDeleteBuffers(1, &obj.mesh.vbo);
        glDeleteBuffers(1, &obj.mesh.ebo);
    }
    CleanupGrid();
    for (auto& pair : gTexCache) glDeleteTextures(1, &pair.second);
    gTexCache.clear();
    withTransform.clear();
    fromNodes.clear();
    glDeleteProgram(shaderProgram);
    glDeleteProgram(gridShaderProgram);
    glfwTerminate();
    return 0;
}