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
    * ================================ VERSION 0.0.3 END ===============================

    * ================================ VERSION 0.0.4 START  ============================
    * fixed dds loading
    * fixed some stuff
    * added more mess
    * json parser obviously (until I dont figure it out, I wont implement it)
    * fixed grid
    * TIRED ASF
    * ================================ VERSION 0.0.4 END ===============================
    
                    
    * ======================= (09/21/2025) V0.0.5 START  ===============================
    * fix (f*ck the version history for this version)
    * significant cleanup as well
    * ================================ VERSION 0.0.5 END ===============================
*/

#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#include "gtc/quaternion.hpp"

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

static const std::string MODELS_ROOT = "C:/drasa_online/work/models/";
static const std::string MESHES_ROOT = "C:/drasa_online/work/meshes/";

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

template<typename T>
T clamp(T v, T lo, T hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}

constexpr double M_PI = 3.14159265358979323846;

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

#define N3_DBG 1
#define N3LOG(x) do{ if(N3_DBG){ std::cout << x << std::endl; } }while(0)

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

struct Options {
    bool ignore_version = true;
    std::string n3filepath;
};

struct Joint {
    int32_t joint_idx = 0;
    int32_t parent_joint_idx = 0;
    std::array<float, 4> pose_translation{ 0,0,0,0 };
    std::array<float, 4> pose_rotation{ 0,0,0,0 };
    std::array<float, 4> pose_scale{ 0,0,0,0 };
    std::string joint_name;
};

struct Node {
    std::string node_name;
    std::string node_type;
    Node* node_parent = nullptr;
    std::vector<Node*> node_children;

    std::array<float, 4> position{ 0.0f,0.0f,0.0f,1.0f };
    std::array<float, 4> rotation{ 0.0f,0.0f,0.0f,0.0f };
    std::array<float, 4> scale{ 1.0f,1.0f,1.0f,1.0f };

    std::string material_name;
    std::string shader_name;

    std::unordered_map<std::string, std::string> shader_textures;
    std::unordered_map<std::string, float>        shader_params1;
    std::unordered_map<std::string, std::array<float, 4>> shader_params4;

    std::unordered_map<int, int> uvSetBySlot;
    std::unordered_map<int, std::array<float, 4>> uvXformBySlot;

    std::string mesh_ressource_id;
    int32_t primitive_group_idx = 0;

    int32_t num_joints = 0;
    std::vector<Joint> joints;

    int32_t num_skin_fragments = 0;
    std::unordered_map<int32_t, std::vector<int32_t>> skin_fragments;

    void SetUvSlot(int slot, int uvSet, float us, float vs, float uo, float vo) {
        uvSetBySlot[slot] = uvSet;
        uvXformBySlot[slot] = { us,vs,uo,vo };
    }
};

struct Reporter {
    enum Type { Info, Warn, Err };
    void report(Type t, const std::string& msg) {
        const char* k = t == Info ? "INFO" : t == Warn ? "WARNING" : "ERROR";
        std::cerr << "[" << k << "] " << msg << "\n";
    }
};

static inline void log_vec4(const char* label, const std::array<float, 4>& v) {
    std::cout << "        " << label << " [" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]\n";
}

class Parser {
public:
    Parser(Reporter& r, const Options& opt) : rep(r), options(opt) {
        canonical_map = {
            {"LDM>", "MODL"}, {">MDL", "MODL"}, {"LDOM","MODL"},
            {"LDM<", "MODL_END"}, {"<MDL","MODL_END"},
            {"DNM>", "MND_START"}, {">MND","MND_START"},
            {"DNM<", "MND_END"}, {"<MND","MND_END"},
            {"_FOE","EOF_"}, {"EOF_","EOF_"},
            {"TRFN","NFRT"}, {"NFRT","NFRT"},
            {"XOBL","LBOX"}, {"LBOX","LBOX"},
            {"IRGP","PGRI"}, {"PGRI","PGRI"},
            {"HSEM","MESH"}, {"SMES","MESH"}, {"MESH","MESH"},
            {"POSI","POSI"}, {"ROTN","ROTN"}, {"SCAL","SCAL"},
            {"MATE","MATE"},
            {"RDHS","SHDR"}, {"SHDR","SHDR"},
            {"STXT","STXT"}, {"TXTS","STXT"},
            {"SFLT","SFLT"}, {"TLFS","SFLT"},
            {"SVEC","SVEC"}, {"SFV4","SVEC"}, {"CEVS","SVEC"},
            {"STUS","STUS"},
            {"SCVU","SCVU"},
            {"NSKF","NSKF"}, {"SFRG","SFRG"},
            {"NJNT","NJNT"}, {"TNJN","NJNT"},
            {"JONT","JONT"}, {"TNOJ","JONT"}
        };
        known_tags = {
            "MODL","MODL_END","MND_START","MND_END","EOF_",
            "MESH","PGRI",
            "POSI","ROTN","SCAL",
            "MATE","SHDR",
            "STXT","SFLT","SVEC",
            "STUS","SCVU",
            "NSKF","SFRG",
            "NJNT","JONT"
        };
    }

    const std::vector<Node*>& getNodes() const { return n3node_list; }

    bool parse_file(const std::string& filepath) {
        this->filepath = filepath;
        std::ifstream f(filepath, std::ios::binary);
        if (!f) { rep.report(Reporter::Err, "Cannot open file"); return false; }
        if (!read_header(f)) return false;

        std::cout << "n3 Version: " << n3version << "\n";

        if (!(n3version == 1 || n3version == 2)) {
            if (options.ignore_version) rep.report(Reporter::Warn, "Unsupported version '" + std::to_string(n3version) + "'");
            else { rep.report(Reporter::Err, "Unsupported version '" + std::to_string(n3version) + "'"); return false; }
        }

        bool done = false; Node* current_node = nullptr; int current_node_idx = -1;

        while (!done) {
            std::string raw; if (!readFourCC(f, raw)) break;
            std::string tag = canonical(raw);

            if (tag == "MODL") {
                std::string t; if (!readFourCC(f, t)) return false;
                n3modeltype = canonical(t);
                n3modelname = readString(f);
                std::cout << "model_type_4cc: '" << n3modeltype << "'\n";
                std::cout << "model_name: '" << n3modelname << "'\n";
                continue;
            }
            if (tag == "MODL_END" || tag == "EOF_") { done = true; break; }

            if (tag == "MND_START") {
                std::string t; if (!readFourCC(f, t)) return false;
                std::string node_type_4cc = canonical(t);
                std::string node_name = readString(f);
                Node* new_node = new Node();
                new_node->node_name = node_name;
                new_node->node_type = node_type_4cc;
                new_node->node_parent = current_node;
                std::cout << "    new_node: " << new_node->node_type << " - " << new_node->node_name << "\n";
                if (current_node) current_node->node_children.push_back(new_node);
                n3node_list.push_back(new_node);
                current_node = new_node;
                current_node_idx += 1;
                continue;
            }
            if (tag == "MND_END") {
                if (current_node) {
                    if (!current_node->material_name.empty())
                        std::cout << "    material: " << current_node->material_name << "\n";
                    if (!current_node->shader_name.empty())
                        std::cout << "    shader: " << current_node->shader_name << "\n";
                    if (!current_node->shader_textures.empty()) {
                        std::cout << "    textures:\n";
                        for (auto& kv : current_node->shader_textures)
                            std::cout << "        " << kv.first << " => " << kv.second << "\n";
                    }
                    if (!current_node->shader_params1.empty()) {
                        std::cout << "    shader_floats:\n";
                        for (auto& kv : current_node->shader_params1)
                            std::cout << "        " << kv.first << " = " << kv.second << "\n";
                    }
                    if (!current_node->shader_params4.empty()) {
                        std::cout << "    shader_vec4:\n";
                        for (auto& kv : current_node->shader_params4) {
                            const auto& v = kv.second;
                            std::cout << "        " << kv.first << " [" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "]\n";
                        }
                    }
                }
                current_node_idx -= 1;
                if (current_node_idx >= 0) current_node = n3node_list[current_node_idx];
                else current_node = nullptr;
                continue;
            }

            if (!current_node) { if (!skip_unknown(f)) break; continue; }

            if (parse_node_tag(f, tag, *current_node)) continue;

            if (!skip_unknown(f)) break;
        }
        return true;
    }

    const std::vector<Node*>& nodes() const { return n3node_list; }
    int version() const { return n3version; }
    const std::string& model_type() const { return n3modeltype; }
    const std::string& model_name() const { return n3modelname; }
    ~Parser() { for (auto* n : n3node_list) delete n; }

private:
    enum Endian { Little, Big } byteOrder = Little;
    Reporter& rep;
    Options options;
    std::string filepath;
    int n3version = 0;
    std::string n3modeltype;
    std::string n3modelname;
    std::vector<Node*> n3node_list;

    std::unordered_map<std::string, std::string> canonical_map;
    std::unordered_set<std::string> known_tags;

    static uint16_t bswap16(uint16_t v) { return uint16_t((v >> 8) | (v << 8)); }
    static uint32_t bswap32(uint32_t v) { return (v >> 24) | ((v >> 8) & 0x0000FF00u) | ((v << 8) & 0x00FF0000u) | (v << 24); }

    bool readU16(std::ifstream& f, uint16_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 2)) return false; if (byteOrder == Big) out = bswap16(out); return true; }
    bool readI32(std::ifstream& f, int32_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 4)) return false; if (byteOrder == Big) out = (int32_t)bswap32((uint32_t)out); return true; }
    bool readU32(std::ifstream& f, uint32_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 4)) return false; if (byteOrder == Big) out = bswap32(out); return true; }
    bool readF32(std::ifstream& f, float& out) { uint32_t u; if (!readU32(f, u)) return false; std::memcpy(&out, &u, 4); return true; }
    bool readI8(std::ifstream& f, int8_t& out) { char c; if (!f.read(&c, 1)) return false; out = (int8_t)c; return true; }
    bool readFourCC(std::ifstream& f, std::string& s) {
        char b[4]; if (!f.read(b, 4)) return false;
        if (byteOrder == Little) s.assign(b, 4); else { char r[4]{ b[3],b[2],b[1],b[0] }; s.assign(r, 4); }
        return true;
    }
    std::string readString(std::ifstream& f) {
        uint16_t n = 0; if (!readU16(f, n)) return {};
        std::string s; s.resize(n);
        if (n > 0) f.read(&s[0], n);
        return s;
    }

    std::string canonical(const std::string& raw) {
        auto it = canonical_map.find(raw);
        if (it != canonical_map.end()) return it->second;
        return raw;
    }

    bool is_known(const std::string& tag) const {
        return known_tags.find(tag) != known_tags.end();
    }

    bool read_header(std::ifstream& f) {
        char magic[4]; if (!f.read(magic, 4)) { rep.report(Reporter::Err, "Failed reading header"); return false; }
        if (std::memcmp(magic, "NEB3", 4) == 0) byteOrder = Little;
        else if (std::memcmp(magic, "3BEN", 4) == 0) byteOrder = Big;
        else { rep.report(Reporter::Err, std::string("Invalid file, unknown fourCC '") + std::string(magic, 4) + "'"); return false; }
        uint32_t ver = 0; if (!readU32(f, ver)) return false;
        if (!(ver == 1 || ver == 2)) {
            uint32_t vs = bswap32(ver);
            if (vs == 1 || vs == 2) { ver = vs; byteOrder = (byteOrder == Little ? Big : Little); }
        }
        n3version = (int)ver;
        return true;
    }

    bool parse_node_tag(std::ifstream& f, const std::string& tag, Node& node) {
        if (tag == "MESH") return parse_mesh(f, node);
        if (tag == "PGRI") return parse_pgri(f, node);
        if (tag == "POSI") return parse_posi(f, node);
        if (tag == "ROTN") return parse_rotn(f, node);
        if (tag == "SCAL") return parse_scal(f, node);
        if (tag == "MATE") return parse_mate(f, node);
        if (tag == "SHDR") return parse_shdr(f, node);
        if (tag == "STXT") return parse_stxt(f, node);
        if (tag == "SFLT") return parse_sflt(f, node);
        if (tag == "SVEC") return parse_svec(f, node);
        if (tag == "STUS") return parse_stus(f, node);
        if (tag == "SCVU") return parse_scvu(f, node);
        if (tag == "NSKF") return parse_nskf(f, node);
        if (tag == "SFRG") return parse_sfrg(f, node);
        if (tag == "NJNT") return parse_njnt(f, node);
        if (tag == "JONT") return parse_jont(f, node);
        return false;
    }

    bool parse_mesh(std::ifstream& f, Node& node) {
        node.mesh_ressource_id = readString(f);
        std::cout << "        MESH " << node.mesh_ressource_id << "\n";
        return true;
    }
    bool parse_pgri(std::ifstream& f, Node& node) {
        int32_t v; if (!readI32(f, v)) return false;
        node.primitive_group_idx = v;
        std::cout << "        PGRI " << node.primitive_group_idx << "\n";
        return true;
    }
    bool parse_posi(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++) if (!readF32(f, node.position[i])) return false;
        log_vec4("POSI", node.position);
        return true;
    }
    bool parse_rotn(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++) if (!readF32(f, node.rotation[i])) return false;
        log_vec4("ROTN", node.rotation);
        return true;
    }
    bool parse_scal(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++) if (!readF32(f, node.scale[i])) return false;
        log_vec4("SCAL", node.scale);
        return true;
    }
    bool parse_mate(std::ifstream& f, Node& node) {
        node.material_name = readString(f);
        std::cout << "        MATE " << node.material_name << "\n";
        return true;
    }
    bool parse_shdr(std::ifstream& f, Node& node) {
        node.shader_name = readString(f);
        std::cout << "        SHDR " << node.shader_name << "\n";
        return true;
    }
    bool parse_stxt(std::ifstream& f, Node& node) {
        std::string tex_type = readString(f);
        std::string tex_name = readString(f);
        node.shader_textures[tex_type] = tex_name;
        std::cout << "        STXT " << tex_type << " => " << tex_name << "\n";
        return true;
    }
    bool parse_sflt(std::ifstream& f, Node& node) {
        std::string pname = readString(f);
        float v; if (!readF32(f, v)) return false;
        node.shader_params1[pname] = v;
        std::cout << "        SFLT " << pname << " = " << v << "\n";
        return true;
    }
    bool parse_svec(std::ifstream& f, Node& node) {
        std::string pname = readString(f);
        std::array<float, 4> v{};
        for (int i = 0; i < 4; i++) if (!readF32(f, v[i])) return false;
        node.shader_params4[pname] = v;
        std::cout << "        SVEC " << pname << " [" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "]\n";
        return true;
    }
    bool parse_stus(std::ifstream& f, Node& node) {
        int32_t pidx; if (!readI32(f, pidx)) return false;
        std::array<float, 4> v{};
        for (int i = 0; i < 4; i++) if (!readF32(f, v[i])) return false;
        node.shader_params4["MLPUVStretch" + std::to_string(pidx)] = v;
        std::cout << "        STUS idx=" << pidx << " [" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << "]\n";
        return true;
    }
    bool parse_scvu(std::ifstream& f, Node& node) {
        int32_t slot; if (!readI32(f, slot)) return false;
        int8_t uvsetRaw; if (!readI8(f, uvsetRaw)) return false;
        (void)readI8(f, uvsetRaw); (void)readI8(f, uvsetRaw); (void)readI8(f, uvsetRaw);
        std::array<float, 4> xform{ 1.0f,1.0f,0.0f,0.0f };
        auto it = node.uvXformBySlot.find(slot);
        if (it != node.uvXformBySlot.end()) xform = it->second;
        node.SetUvSlot(slot, (int)(uint8_t)uvsetRaw, xform[0], xform[1], xform[2], xform[3]);
        std::cout << "        SCVU slot=" << slot << " uvSet=" << (int)(uint8_t)uvsetRaw << "\n";
        return true;
    }
    bool parse_nskf(std::ifstream& f, Node& node) {
        int32_t v; if (!readI32(f, v)) return false;
        node.num_skin_fragments = v;
        std::cout << "        NSKF num_skin_fragments=" << v << "\n";
        return true;
    }
    bool parse_sfrg(std::ifstream& f, Node& node) {
        int32_t group_idx; if (!readI32(f, group_idx)) return false;
        int32_t num_joints; if (!readI32(f, num_joints)) return false;
        std::vector<int32_t> joints(num_joints);
        for (int i = 0; i < num_joints; i++) if (!readI32(f, joints[i])) return false;
        node.skin_fragments[group_idx] = std::move(joints);
        std::cout << "        SFRG group=" << group_idx << " joints=" << num_joints << "\n";
        return true;
    }
    bool parse_njnt(std::ifstream& f, Node& node) {
        int32_t v; if (!readI32(f, v)) return false;
        node.num_joints = v;
        std::cout << "        NJNT num_joints=" << v << "\n";
        return true;
    }
    bool parse_jont(std::ifstream& f, Node& node) {
        Joint j;
        if (!readI32(f, j.joint_idx)) return false;
        if (!readI32(f, j.parent_joint_idx)) return false;
        for (int i = 0; i < 4; i++) if (!readF32(f, j.pose_translation[i])) return false;
        for (int i = 0; i < 4; i++) if (!readF32(f, j.pose_rotation[i])) return false;
        for (int i = 0; i < 4; i++) if (!readF32(f, j.pose_scale[i])) return false;
        j.joint_name = readString(f);
        node.joints.push_back(j);
        std::cout << "        JONT idx=" << j.joint_idx << " parent=" << j.parent_joint_idx << " name=" << j.joint_name << "\n";
        log_vec4("            pose_translation", j.pose_translation);
        log_vec4("            pose_rotation", j.pose_rotation);
        log_vec4("            pose_scale", j.pose_scale);
        return true;
    }

    bool skip_unknown(std::ifstream& f) {
        std::streampos base = f.tellg();
        std::streampos endpos;
        { auto save = f.tellg(); f.seekg(0, std::ios::end); endpos = f.tellg(); f.seekg(save); }
        const int64_t MAX_SCAN = 4 * 1024 * 1024;

        for (int64_t i = 0; i < MAX_SCAN && (base + std::streamoff(i)) < endpos; ++i) {
            f.seekg(base + std::streamoff(i));
            std::string t;
            if (!readFourCC(f, t)) return false;
            std::string ct = canonical(t);
            if (is_known(ct)) {
                f.seekg(base + std::streamoff(i));
                return true;
            }
        }
        return false;
    }
};

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
static inline float ub_to_n11(uint8_t b) { return (float(b) * (1.0f / 255.0f)) * 2.0f - 1.0f; }
static inline float ub_to_u01(uint8_t b) { return float(b) * (1.0f / 255.0f); }
static inline float s_to_n11(int16_t s) { return s < 0 ? float(s) / 32768.0f : float(s) / 32767.0f; }
static inline float s_to_u01(int16_t s) { return (float(s) + 32768.0f) / 65535.0f; }
static inline float s2_to_uv(int16_t s) { return float(s) / 8191.0f; }
static inline float snorm_to_u01(int16_t s) { return 0.5f * (s_to_n11(s) + 1.0f); }
static inline float fps2float(int16_t s) { return float(s) / 8191.0f; }

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

bool LoadNVX2(const std::string& path, Mesh& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { errout << "[NVX2] cannot open " << path << " "; return false; }

    char magic[4]; f.read(magic, 4); if (!f)return false;
    uint32_t numGroups, numVertices, vertexWidthFloats, numTrianglesOrIndices, numEdges, compMask;
    f.read((char*)&numGroups, 4);
    f.read((char*)&numVertices, 4);
    f.read((char*)&vertexWidthFloats, 4);
    f.read((char*)&numTrianglesOrIndices, 4);
    f.read((char*)&numEdges, 4);
    f.read((char*)&compMask, 4);

    out.groups.assign(numGroups, {});
    if (numGroups) { f.read((char*)out.groups.data(), numGroups * sizeof(Nvx2Group)); if (!f)return false; }

    const uint32_t stride = vertexWidthFloats * 4;
    std::unordered_map<uint32_t, int> offsets;
    int off = 0;
    auto add = [&](uint32_t c, int sz) { if (compMask & c) { offsets[c] = off; off += sz; } };
    add(Coord, 12);
    add(Normal, 12); add(NormalUB4N, 4);
    add(Uv0, 8); add(Uv0S2, 4);
    add(Uv1, 8); add(Uv1S2, 4);
    add(Uv2, 8); add(Uv2S2, 4);
    add(Uv3, 8); add(Uv3S2, 4);
    add(Color, 16); add(ColorUB4N, 4);
    add(Tangent, 12); add(TangentUB4N, 4);
    add(Binormal, 12); add(BinormalUB4N, 4);
    add(Weights, 16); add(WeightsUB4N, 4);
    add(JIndices, 16); add(JIndicesUB4, 4);

    std::vector<uint8_t> vbuf(size_t(numVertices) * stride);
    if (!vbuf.empty()) { f.read((char*)vbuf.data(), vbuf.size()); if (!f)return false; }

    out.verts.assign(numVertices, {});
    for (uint32_t i = 0; i < numVertices; i++) {
        const uint8_t* base = vbuf.data() + size_t(i) * stride;
        ObjVertex v{}; v.nz = 1; v.cr = v.cg = v.cb = v.ca = 1; v.w0 = 1;
        if (offsets.count(Coord)) {
            const float* p = (const float*)(base + offsets[Coord]);
            v.px = p[0]; v.py = p[1]; v.pz = p[2];
        }
        if (offsets.count(Normal)) {
            const float* n = (const float*)(base + offsets[Normal]);
            v.nx = n[0]; v.ny = n[1]; v.nz = n[2];
        }
        else if (offsets.count(NormalUB4N)) {
            const uint8_t* n = base + offsets[NormalUB4N];
            v.nx = ub_to_n11(n[0]); v.ny = ub_to_n11(n[1]); v.nz = ub_to_n11(n[2]);
            float len = std::sqrt(v.nx * v.nx + v.ny * v.ny + v.nz * v.nz); if (len > 0) { float inv = 1.0f / len; v.nx *= inv; v.ny *= inv; v.nz *= inv; }
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
        if (offsets.count(Uv0)) { const float* u = (const float*)(base + offsets[Uv0]); v.u0 = u[0]; v.v0 = u[1]; }
        else if (offsets.count(Uv0S2)) { const int16_t* u = (const int16_t*)(base + offsets[Uv0S2]); v.u0 = fps2float(u[0]); v.v0 = fps2float(u[1]); }
        if (offsets.count(Uv1)) { const float* u = (const float*)(base + offsets[Uv1]); v.u1 = u[0]; v.v1 = u[1]; }
        else if (offsets.count(Uv1S2)) { const int16_t* u = (const int16_t*)(base + offsets[Uv1S2]); v.u1 = fps2float(u[0]); v.v1 = fps2float(u[1]); }
        if (offsets.count(Color)) { const float* c = (const float*)(base + offsets[Color]); v.cr = c[0]; v.cg = c[1]; v.cb = c[2]; v.ca = c[3]; }
        else if (offsets.count(ColorUB4N)) { const uint8_t* c = base + offsets[ColorUB4N]; v.cr = c[0] / 255.0f; v.cg = c[1] / 255.0f; v.cb = c[2] / 255.0f; v.ca = c[3] / 255.0f; }
        if (offsets.count(Weights)) { const float* w = (const float*)(base + offsets[Weights]); v.w0 = w[0]; v.w1 = w[1]; v.w2 = w[2]; v.w3 = w[3]; }
        else if (offsets.count(WeightsUB4N)) {
            const uint8_t* w = base + offsets[WeightsUB4N];
            v.w0 = w[0] / 255.0f; v.w1 = w[1] / 255.0f; v.w2 = w[2] / 255.0f; v.w3 = w[3] / 255.0f;
            float s = v.w0 + v.w1 + v.w2 + v.w3; if (s > 0) { float inv = 1.0f / s; v.w0 *= inv; v.w1 *= inv; v.w2 *= inv; v.w3 *= inv; }
        }
        if (offsets.count(JIndices)) { const float* jf = (const float*)(base + offsets[JIndices]); v.j0 = (uint8_t)jf[0]; v.j1 = (uint8_t)jf[1]; v.j2 = (uint8_t)jf[2]; v.j3 = (uint8_t)jf[3]; }
        else if (offsets.count(JIndicesUB4)) { const uint8_t* jb = base + offsets[JIndicesUB4]; v.j0 = jb[0]; v.j1 = jb[1]; v.j2 = jb[2]; v.j3 = jb[3]; }
        out.verts[i] = v;
    }

    uint32_t totalTris = 0; for (auto& g : out.groups) totalTris += g.numTriangles;
    uint32_t expectedIndexCount = totalTris ? totalTris * 3 : numTrianglesOrIndices;

    out.idx.clear();
    if (expectedIndexCount) {
        std::streampos pos = f.tellg(); f.seekg(0, std::ios::end);
        size_t remain = size_t(f.tellg() - pos); f.seekg(pos);
        size_t need16 = size_t(expectedIndexCount) * 2, need32 = size_t(expectedIndexCount) * 4;
        if ((remain >= need32 && remain % 4 == 0 && numVertices >= 65536) || remain == need32) {
            std::vector<uint32_t> tmp(expectedIndexCount);
            f.read((char*)tmp.data(), need32); if (!f)return false;
            out.idx = std::move(tmp);
        }
        else {
            std::vector<uint16_t> tmp(expectedIndexCount);
            f.read((char*)tmp.data(), need16); if (!f)return false;
            out.idx.resize(expectedIndexCount);
            for (uint32_t i = 0; i < expectedIndexCount; i++) out.idx[i] = tmp[i];
        }
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

uniform int HasBump;
uniform int HasSpec;

void main() {
    vec4 diffTex = texture(DiffMap, vUV0);
    if (diffTex.a < 0.5) discard;

    vec3 N = normalize(vNormal);
    if (HasBump == 1) {
        vec3 mapN = texture(BumpMap, vUV0).xyz * 2.0 - 1.0;
        N = normalize(mix(N, mapN, 1.0));
    }

    vec3 V = normalize(CameraPos - vPos);
    vec3 L = normalize(-SunDir);
    vec3 H = normalize(L + V);

    float ndl = max(dot(N, L), 0.0);
    float ndh = max(dot(N, H), 0.0);

    vec3 albedo = diffTex.rgb;
    vec3 specC = (HasSpec == 1) ? texture(SpecMap, vUV1).rgb : vec3(0.0);
    float spec = (ndl > 0.0) ? pow(ndh, SpecularPower) * SpecularIntensity : 0.0;
    vec3 emissive = texture(EmissiveMap, vUV1).rgb * EmissiveIntensity;

    vec3 ambient = AmbientColor * albedo;
    vec3 diffuse = SunColor * SunIntensity * ndl * albedo;
    vec3 specular = SunColor * SunIntensity * spec * specC;

    vec3 color = ambient + diffuse + specular + emissive;
    FragColor = vec4(color, diffTex.a);
}
)";


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

class Camera {
    enum class CameraMovement {
        FORWARD,
        BACKWARD,
        LEFT,
        RIGHT,
        UP,
        DOWN
    };
public:
    Camera(std::string  name, glm::vec3 position, glm::vec3 forward, glm::vec3 up,
        float yaw, float pitch, float moveSpeed,
        float mouseSensitivity, float fov, float nearPlane, float farPlane)
        : m_position(position),
        m_forwardVec(forward),
        m_upVec(up),
        name(std::move(name)),
        m_yaw(yaw),
        m_pitch(pitch),
        m_movementSpeed(moveSpeed),
        m_mouseSensitivity(mouseSensitivity),
        m_fov(fov),
        m_nearPlane(nearPlane),
        m_farPlane(farPlane)
    {
        m_yaw = glm::degrees(atan2(m_forwardVec.z, m_forwardVec.x));
        m_pitch = glm::degrees(asin(m_forwardVec.y));

        updateProjectionMatrix();
        updateVectors();
        updateViewMatrix();
    }
    glm::mat4 getViewMatrix() {
        return m_viewMatrix;
    }

    void setProjectionMatrix(const glm::mat4& projection) {
        m_ProjectionMatrix = projection;
    }

    const glm::mat4& getProjectionMatrix() const {
        return m_ProjectionMatrix;
    }

    std::string getName() const {
        return name;
    }

    void processKeyboard(CameraMovement direction, float deltaTime) {
        float velocity = m_movementSpeed * deltaTime;
        if (direction == CameraMovement::FORWARD)  m_position += m_forwardVec * velocity;
        if (direction == CameraMovement::BACKWARD) m_position -= m_forwardVec * velocity;
        if (direction == CameraMovement::LEFT)     m_position -= glm::normalize(glm::cross(m_forwardVec, m_upVec)) * velocity;
        if (direction == CameraMovement::RIGHT)    m_position += glm::normalize(glm::cross(m_forwardVec, m_upVec)) * velocity;
        if (direction == CameraMovement::UP)       m_position.y += velocity;
        if (direction == CameraMovement::DOWN)     m_position.y -= velocity;
        updateViewMatrix();
    }

    void processMouseMovement(float xoffset, float yoffset, GLboolean constrainPitch = true) {
        xoffset *= m_mouseSensitivity;
        yoffset *= m_mouseSensitivity;

        m_yaw += xoffset;
        m_pitch += yoffset;

        if (constrainPitch) {
            if (m_pitch > 89.0f)  m_pitch = 89.0f;
            if (m_pitch < -89.0f) m_pitch = -89.0f;
        }
        updateVectors();
        updateViewMatrix();
    }

    void processMouseScroll(float yoffset) {
        m_fov -= yoffset;
        if (m_fov < 1.0f)  m_fov = 1.0f;
        if (m_fov > 45.0f) m_fov = 45.0f;
        updateProjectionMatrix();
    }

    void handleInput(GLFWwindow* window, float deltaTime) {
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) processKeyboard(CameraMovement::FORWARD, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) processKeyboard(CameraMovement::BACKWARD, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) processKeyboard(CameraMovement::LEFT, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) processKeyboard(CameraMovement::RIGHT, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) processKeyboard(CameraMovement::UP, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS) processKeyboard(CameraMovement::DOWN, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) processKeyboard(CameraMovement::UP, deltaTime);
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) processKeyboard(CameraMovement::DOWN, deltaTime);
    }

    static void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
        Camera* cam = reinterpret_cast<Camera*>(glfwGetWindowUserPointer(window));
        if (!cam) return;
        static float lastX = 0.0f;
        static float lastY = 0.0f;
        static bool firstMouse = true;
        if (firstMouse) {
            lastX = (float)xpos;
            lastY = (float)ypos;
            firstMouse = false;
        }
        float xoffset = (float)xpos - lastX;
        float yoffset = lastY - (float)ypos;
        lastX = (float)xpos;
        lastY = (float)ypos;
        cam->processMouseMovement(xoffset, yoffset);
    }

    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
        Camera* cam = reinterpret_cast<Camera*>(glfwGetWindowUserPointer(window));
        if (!cam) return;
        cam->processMouseScroll((float)yoffset);
    }

    void processScroll(double yOffset) {
        float scrollSpeed = 1.5f;
        m_position += scrollSpeed * static_cast<float>(yOffset) * m_forwardVec;
        updateViewMatrix();
    }

    float getFov() const {
        return m_fov;
    }

    glm::vec3 getPosition() const {
        return m_position;
    }

    float getYaw() const {
        return m_yaw;
    }

    float getPitch() const {
        return m_pitch;
    }

    float setYaw(float newYaw) {
        return m_yaw = newYaw;
    }

    float setPitch(float newPitch) {
        return m_pitch = newPitch;
    }

    float getFOV() const {
        return m_fov;
    }

    void setFOV(float fov) {
        m_fov = fov;
        updateProjectionMatrix();
    }

    float getNearPlane() const {
        return m_nearPlane;
    }

    void setNearPlane(float nearPlane) {
        m_nearPlane = nearPlane;
        updateProjectionMatrix();
    }

    float getFarPlane() const {
        return m_farPlane;
    }

    void setFarPlane(float farPlane) {
        m_farPlane = farPlane;
        updateProjectionMatrix();
    }

    void setName(const std::string& newName) {
        name = newName;
    }

    void setPosition(const glm::vec3& position) {
        m_position = position;
        updateViewMatrix();
    }

    glm::vec3 getForwardVector() const {
        return m_forwardVec;
    }

    glm::vec3 getUpVector() const {
        return m_upVec;
    }

    void setUpVector(const glm::vec3& upVec) {
        m_upVec = upVec;
        updateViewMatrix();
    }

private:
    void updateVectors() {
        glm::vec3 forward;
        forward.x = cos(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
        forward.y = sin(glm::radians(m_pitch));
        forward.z = sin(glm::radians(m_yaw)) * cos(glm::radians(m_pitch));
        m_forwardVec = glm::normalize(forward);

        right = glm::normalize(glm::cross(m_forwardVec, glm::vec3(0.0f, 1.0f, 0.0f)));
        m_upVec = glm::normalize(glm::cross(right, m_forwardVec));
    }
    void updateViewMatrix() {
        glm::vec3 target = m_position + m_forwardVec;
        m_viewMatrix = glm::lookAt(m_position, target, m_upVec);
    }

    void updateProjectionMatrix() {
        m_ProjectionMatrix = glm::perspective(glm::radians(m_fov), 16.0f / 9.0f, m_nearPlane, m_farPlane);
    }

    std::string name;
    glm::mat4 m_ProjectionMatrix{};
    glm::mat4 m_viewMatrix{};
    glm::vec3 right{};
    glm::vec3 m_position;
    glm::vec3 m_forwardVec;
    glm::vec3 m_upVec;

    float m_yaw;
    float m_pitch;
    float m_movementSpeed;
    float m_mouseSensitivity;
    float m_fov;
    float m_nearPlane;
    float m_farPlane;
};


std::vector<DrawCmd> BuildDrawsWithTransform(const Parser& parser, const glm::vec3& pos, const glm::vec3& rot, const glm::vec3& scale) {
    std::vector<DrawCmd> outv;
    glm::mat4 trs = composeTRS(pos, rot, scale);
    for (auto* node : parser.getNodes()) {
        if (node->mesh_ressource_id.empty()) continue;
        DrawCmd dc;
        dc.group = node->primitive_group_idx;
        glm::mat4 nodeLocal;
        buildTransformMatrix(&nodeLocal[0][0], Vec3{ node->position[0],node->position[1],node->position[2] }, Vec4{ node->rotation[0],node->rotation[1],node->rotation[2],node->rotation[3] }, Vec3{ node->scale[0],node->scale[1],node->scale[2] });
        glm::mat4 finalM = trs * nodeLocal; matCopy(dc.worldMatrix, finalM);
        dc.tex[0] = loadNebTex(node->shader_textures, "DiffMap0");
        dc.tex[1] = loadNebTex(node->shader_textures, "SpecMap0");
        dc.tex[2] = loadNebTex(node->shader_textures, "BumpMap0");
        dc.tex[3] = loadNebTex(node->shader_textures, "EmissiveMap0");
        for (int i = 0; i < 4; i++) { dc.has[i] = (dc.tex[i] != 0); if (!dc.has[i]) dc.tex[i] = (i == 0) ? gWhiteTex : gBlackTex; }
        std::string meshPath = node->mesh_ressource_id;
        if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
        meshPath = MESHES_ROOT + meshPath;
        if (!LoadNVX2(meshPath, dc.mesh)) continue;
        dc.name = meshPath;
        for (int s = 0; s < 4; s++) { dc.uvXform[s][0] = 1; dc.uvXform[s][1] = 1; dc.uvXform[s][2] = 0; dc.uvXform[s][3] = 0; dc.uvSet[s] = 0; }
        outv.push_back(std::move(dc));
    }
    return outv;
}


static inline std::string strlower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return (char)std::tolower(c); });
    return s;
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

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) return -1;
    glViewport(0, 0, (int)sizeX, (int)sizeY);

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

    Camera camera(
        "Main Camera",
        glm::vec3(0.0f, 0.0f, 5.0f),
        glm::vec3(0.0f, 0.0f, -1.0f),
        glm::vec3(0.0f, 1.0f, 0.0f),
        -90.0f, 0.0f,
        5.0f,
        0.1f,
        45.0f,
        0.1f, 1500.0f
    );

    float lastFrame = 0.0f;

    glfwSetWindowUserPointer(window, &camera);
    glfwSetCursorPosCallback(window, Camera::mouse_callback);
    glfwSetScrollCallback(window, Camera::scroll_callback);
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    auto appendN3WTransform = [&](const std::string& path,
        const glm::vec3& pos,
        const glm::vec3& rot,
        const glm::vec3& scale)
        {
            Reporter rep;
            Options opt;
            Parser parser(rep, opt);

            if (!parser.parse_file(path)) return;

            auto draws = BuildDrawsWithTransform(parser, pos, rot, scale);
            for (auto& dc : draws) withTransform.push_back(std::move(dc));
        };

    appendN3WTransform("C:/drasa_online/work/models/t001_hub/arch_house_01.n3",
        { 0.0,0.0,0.0 }, { 0,0,0 }, { 1,1,1 });
    appendN3WTransform("C:/drasa_online/work/models/t001_hub/arch_house_roof_01.n3",
        { 0.0,0.0,0.0 }, { 0,0,0 }, { 1,1,1 }); 

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

    glUseProgram(0);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        float deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        camera.handleInput(window, deltaTime);

        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);

        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(camera.getProjectionMatrix()));
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(camera.getViewMatrix()));

        glm::vec3 camPos = camera.getPosition();
        glUniform3fv(glGetUniformLocation(shaderProgram, "CameraPos"), 1, glm::value_ptr(camPos));
        glUniform3fv(glGetUniformLocation(shaderProgram, "PointPos"), 1, glm::value_ptr(camPos));


        GLint locUseUV1 = glGetUniformLocation(shaderProgram, "UseUV1");
        GLint locDiffMap = glGetUniformLocation(shaderProgram, "DiffMap");

        for (auto& obj : withTransform) {
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, obj.worldMatrix);

            glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, obj.tex[0]);
            glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, obj.tex[1]);
            glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, obj.tex[2]);
            glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, obj.tex[3]);

            glUniform1i(glGetUniformLocation(shaderProgram, "HasSpec"), obj.has[1] ? 1 : 0);
            glUniform1i(glGetUniformLocation(shaderProgram, "HasBump"), obj.has[2] ? 1 : 0);

            obj.applyUVTransforms(shaderProgram);

            glBindVertexArray(obj.mesh.vao);

            if (!obj.mesh.idx.empty()) {
                if (obj.group >= 0 && obj.group < (int)obj.mesh.groups.size()) {
                    const Nvx2Group& g = obj.mesh.groups[obj.group];
                    const GLsizei count = (GLsizei)g.indexCount();
                    const GLsizei offset = (GLsizei)(g.firstIndex() * sizeof(uint32_t));
                    glDrawElements(GL_TRIANGLES, count, GL_UNSIGNED_INT, (void*)(intptr_t)offset);
                }
                else {
                    glDrawElements(GL_TRIANGLES, (GLsizei)obj.mesh.idx.size(), GL_UNSIGNED_INT, 0);
                }
            }

            glBindVertexArray(0);
        }

        glUseProgram(0);
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    for (auto& obj : withTransform) {
        glDeleteVertexArrays(1, &obj.mesh.vao);
        glDeleteBuffers(1, &obj.mesh.vbo);
        glDeleteBuffers(1, &obj.mesh.ebo);
    }

    for (auto& pair : gTexCache) glDeleteTextures(1, &pair.second);
    gTexCache.clear();
    withTransform.clear();
    glDeleteProgram(shaderProgram);
    glfwTerminate();
    return 0;
}