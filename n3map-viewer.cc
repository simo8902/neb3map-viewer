//PRODBYSIMO 101625
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#undef APIENTRY
#define APIENTRY __stdcall
#include <wincodec.h>
#include <combaseapi.h>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <DirectXTex/DirectXTex.h>

#define GLM_ENABLE_EXPERIMENTAL
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"
#include "gtc/type_ptr.hpp"
#include "gtc/quaternion.hpp"
#include "gtx/quaternion.hpp"
#include "gtx/string_cast.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <array>
#include <random>
#include <unordered_set>
#include <sstream>
#include <set>
#include <filesystem>
#include <fstream>
#include <string>
#include <cctype>

static const std::string& MODELS_ROOT = "C:/drasa_online/work/models/";
static const std::string& MESHES_ROOT = "C:/drasa_online/work/meshes/";
static const std::string& ANIMS_ROOT = "C:/drasa_online/work/anims/"; // stripped
static const std::string& MAP_ROOT = "C:/drasa_online/work/maps/"; // stripped

void InitDeferred(int w, int h);
void ResizeDeferred(int w, int h);

int gFBWidth = 1200, gFBHeight = 800;
GLuint gBuffer = 0, gPosition = 0, gNormal = 0, gAlbedoSpec = 0, gDepth = 0;
inline GLFWwindow* window = nullptr;
inline GLuint shaderProgram = 0;
GLuint gSamplerRepeat = 0;
GLuint gSamplerClamp = 0;

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

inline GLuint quadVAO = 0;
inline GLuint quadVBO = 0;

constexpr double M_PI = 3.14159265358979323846;

struct ObjVertex { float px, py, pz; float nx, ny, nz; float tx, ty, tz; float bx, by, bz; float u0, v0; float u1, v1; float cr, cg, cb, ca; float w0, w1, w2, w3; uint8_t j0, j1, j2, j3; };

#pragma pack(push, 1)
struct Nvx2Group {
    uint32_t firstVertex;
    uint32_t numVertices;
    uint32_t firstTriangle;
    uint32_t numTriangles;
    uint32_t firstEdge;
    uint32_t numEdges;

    uint32_t firstIndex()   const { return firstTriangle * 3; }
    uint32_t indexCount()   const { return numTriangles * 3; }
    uint32_t baseVertex()   const { return firstVertex; }
};
#pragma pack(pop)

struct Mesh {
    std::vector<ObjVertex> verts;
    std::vector<uint32_t> idx;
    std::vector<Nvx2Group> groups;
    GLuint vao = 0, vbo = 0, ebo = 0;
};

struct DrawCmd {
    Mesh mesh;
    std::string name;
    std::string nodeName;
    std::string mat;
    std::string shdr;
    std::unordered_map<std::string, std::string> shader_textures;

    GLuint tex[8];
    bool   has[8];
    glm::mat4 worldMatrix{ 1.0f };
    glm::mat4 rootMatrix{ 1.0f };
    int    group = -1;
    float textureTile = 1.0f;
};

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

struct Options {
    bool ignore_version = true;
    std::string n3filepath;
};

struct Node {
    std::string node_name;
    std::string node_type;
    Node* node_parent = nullptr;
    std::vector<Node*> node_children;

    std::array<float, 4> position{ 0.0f, 0.0f, 0.0f, 1.0f };
    std::array<float, 4> rotation{ 0.0f, 0.0f, 0.0f, 1.0f };
    std::array<float, 4> scale{ 1.0f, 1.0f, 1.0f, 1.0f };

    std::unordered_map<std::string, std::string> shader_textures;
    std::string mesh_ressource_id;
    int32_t primitive_group_idx = 0;
    std::string mat;
    std::string shdr;

};

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

    const glm::mat4& getProjectionMatrix() const {
        return m_ProjectionMatrix;
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


    glm::vec3 getPosition() const {
        return m_position;
    }

    void updateAspectRatio(int width, int height) {
        if (height > 0) {
            m_aspectRatio = static_cast<float>(width) / static_cast<float>(height);
            updateProjectionMatrix();
        }
    }

    float getAspectRatio() const {
        return m_aspectRatio;
    }

    void setPerspective(float fovDeg, float aspect, float nearPlane, float farPlane) {
        m_fov = fovDeg;
        m_aspectRatio = aspect;
        m_nearPlane = nearPlane;
        m_farPlane = farPlane;
        updateProjectionMatrix();
    }

    float getFov() const { return m_fov; }
    void setFov(float fovDeg) { m_fov = fovDeg; updateProjectionMatrix(); }

    float getNearPlane() const { return m_nearPlane; }
    float getFarPlane()  const { return m_farPlane; }
    void setNearPlane(float v) { m_nearPlane = v; updateProjectionMatrix(); }
    void setFarPlane(float v) { m_farPlane = v; updateProjectionMatrix(); }

    void setAspectRatio(float aspect) { m_aspectRatio = aspect; updateProjectionMatrix(); }
    glm::vec4 getPerspective() const { return glm::vec4(m_fov, m_aspectRatio, m_nearPlane, m_farPlane); }

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
        m_ProjectionMatrix = glm::perspective(glm::radians(m_fov), m_aspectRatio, m_nearPlane, m_farPlane);
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
    float m_aspectRatio = 16.0f / 9.0f;
};


Camera camera = Camera(
    "MainCamera",
    glm::vec3(0.0f, 0.0f, 3.0f),
    glm::vec3(0.0f, 0.0f, -1.0f),
    glm::vec3(0.0f, 1.0f, 0.0f),
    0.0f, 0.0f, 15.0f, 0.1f, 45.0f, 0.1f, 1000.0f
);

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    gFBWidth = width;
    gFBHeight = height;

    glViewport(0, 0, gFBWidth, gFBHeight);

    Camera* cam = static_cast<Camera*>(glfwGetWindowUserPointer(window));
    if (cam) cam->updateAspectRatio(gFBWidth, gFBHeight);

    ResizeDeferred(gFBWidth, gFBHeight);
}

inline void InitGLFW() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SRGB_CAPABLE, GL_TRUE);

    int sizeX = 1200, sizeY = 800;

    // STRIPPED! for full pay me
    window = glfwCreateWindow(sizeX, sizeY, "SIMO'S NEB3 VIEWER", nullptr, nullptr);
    if (!window) { glfwTerminate(); throw std::runtime_error("GLFW window failed"); }

    glfwMakeContextCurrent(window);
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        throw std::runtime_error("Failed to init GLAD");

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glEnable(GL_FRAMEBUFFER_SRGB);
    glDisable(GL_BLEND);

    glfwGetFramebufferSize(window, &gFBWidth, &gFBHeight);
    glViewport(0, 0, gFBWidth, gFBHeight);

    InitDeferred(gFBWidth, gFBHeight);

    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
}

static inline void log_vec4(const char* label, const std::array<float, 4>& v) {
    std::cout << "        " << label << " [" << v[0] << ", " << v[1] << ", " << v[2] << ", " << v[3] << "]\n";
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
    glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, u1));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, tx));

    glEnableVertexAttribArray(7);
    glVertexAttribPointer(7, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, bx));


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

GLuint shaderStandard;      // shd:standard
GLuint shaderLighting;
GLuint shaderEnvironment = 0;
GLuint shaderSimpleLayer;
GLuint shaderDecal;

void InitSHDR() {
    auto make = [&](const char* vs, const char* fs) -> GLuint {
        GLuint v = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(v, 1, &vs, nullptr);
        glCompileShader(v);
        CheckShader(v, "VS");

        GLuint f = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(f, 1, &fs, nullptr);
        glCompileShader(f);
        CheckShader(f, "FS");

        GLuint prog = glCreateProgram();
        glAttachShader(prog, v);
        glAttachShader(prog, f);
        glLinkProgram(prog);
        CheckProgram(prog);

        glDeleteShader(v);
        glDeleteShader(f);
        return prog;
        };

    const char* deferred_vs = R"(
#version 330 core
layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 texcoord0;
layout(location = 3) in vec4 aJointWeights;
layout(location = 4) in ivec4 aJointIndices;
layout(location = 5) in vec2 texcoord1;
layout(location = 6) in vec3 tangent;
layout(location = 7) in vec3 binormal;

uniform mat4 projection, view, model;
uniform mat4 JointMatrices[128];
uniform int UseSkinning;

out vec3 vWorldPos;
out vec3 vNormalWS;
out vec2 vUV;
out vec2 vUV1;
out vec3 vTangentWS;
out vec3 vBinormalWS;

void main() {
vec4 localPos = vec4(position, 1.0);
vec3 localNorm = normal;
vec3 localTan = tangent;
vec3 localBin = binormal;
if (UseSkinning > 0) {
    vec4 w = aJointWeights;
    float s = w.x + w.y + w.z + w.w;
    if (s > 0.0) w /= s;
    mat4 skinMat =
        w.x * JointMatrices[aJointIndices.x] +
        w.y * JointMatrices[aJointIndices.y] +
        w.z * JointMatrices[aJointIndices.z] +
        w.w * JointMatrices[aJointIndices.w];
    localPos = skinMat * localPos;
    localNorm = mat3(skinMat) * localNorm;
    localTan  = mat3(skinMat) * localTan;
    localBin  = mat3(skinMat) * localBin;
}

vec4 wpos = model * localPos;
vWorldPos = wpos.xyz;

mat3 nM = mat3(transpose(inverse(model)));
vNormalWS   = normalize(nM * localNorm);
vTangentWS  = normalize(nM * localTan);
vBinormalWS = normalize(nM * localBin);

vUV  = texcoord0;
vUV1 = texcoord1;

gl_Position = projection * (view * wpos);
}

)";

    const char* deferred_fs = R"(
#version 330 core

uniform sampler2D DiffMap;
uniform sampler2D SpecMap;
uniform sampler2D BumpMap;
uniform sampler2D EmsvMap; 

uniform int   alphaTest;      
uniform float alphaCutoff; 
uniform int   twoSided;      
uniform int   isFlatNormal;

in vec3 vWorldPos;
in vec3 vNormalWS;
in vec2 vUV;
in vec3 vTangentWS;
in vec3 vBinormalWS;

layout(location=0) out vec4 gPosition;
layout(location=1) out vec4 gNormal;
layout(location=2) out vec4 gAlbedoSpec;

vec3 SampleNormalDXT5(vec4 s){
    vec2 xy = s.wy * 2.0 - 1.0;
    float z = sqrt(max(0.0, 1.0 - dot(xy,xy)));
    return vec3(xy,z);
}

void main(){
    vec4 d = texture(DiffMap, vUV);
    if (alphaTest == 1 && d.a < alphaCutoff) discard;

    float sR = texture(SpecMap, vUV).r;
    vec4  b  = texture(BumpMap, vUV);

    vec3 N = normalize(vNormalWS);
    vec3 T = normalize(vTangentWS);
    vec3 B = normalize(vBinormalWS);

    if (twoSided == 1 && !gl_FrontFacing) {
        N = -N; B = -B;
    }

    mat3 TBN = mat3(T,B,N);
    vec3 tn  = (isFlatNormal == 1) ? vec3(0,0,1) : SampleNormalDXT5(b);
    vec3 Nws = normalize(TBN * tn);

    gPosition   = vec4(vWorldPos, 1.0);
    gNormal     = vec4(Nws * 0.5 + 0.5, 1.0);
    gAlbedoSpec = vec4(d.rgb, sR);
}
)";

    const char* lighting_vs = R"(
#version 330 core
layout(location=0) in vec2 aPos;
layout(location=1) in vec2 aUV;
out vec2 TexCoords;
void main(){
    TexCoords = aUV;
    gl_Position = vec4(aPos,0.0,1.0);
}
)";

    const char* lighting_fs = R"(
#version 330 core
in vec2 TexCoords;
out vec4 FragColor;
uniform sampler2D gAlbedoSpec;
void main(){
    vec4 col = texture(gAlbedoSpec, TexCoords);
    FragColor = length(col.rgb) < 0.001 ? vec4(1,0,1,1) : vec4(col.rgb,1);
}
)";

    const char* simplelayer_gbuffer_fs = R"(
#version 330 core
uniform sampler2D DiffMap0;
uniform sampler2D SpecMap0;
uniform sampler2D BumpMap0;
uniform sampler2D DiffMap2;
uniform sampler2D SpecMap1;
uniform sampler2D BumpMap1;
uniform sampler2D DiffMap3;

uniform int twoSided;
uniform int isFlat0;
uniform int isFlat1;

in vec3 vWorldPos;
in vec3 vNormalWS;
in vec2 vUV;
in vec3 vTangentWS;
in vec3 vBinormalWS;

layout(location=0) out vec4 gPosition;
layout(location=1) out vec4 gNormal;
layout(location=2) out vec4 gAlbedoSpec;

vec3 SampleNormalDXT5(vec4 s){
    vec2 xy = s.wy * 2.0 - 1.0;
    float z = sqrt(max(0.0, 1.0 - dot(xy, xy)));
    return vec3(xy, z);
}

vec3 BlendNormals(vec3 n0, vec3 n1, float w){
    n0.z = max(n0.z, 1e-3);
    n1.z = max(n1.z, 1e-3);
    return normalize(mix(n0, n1, w));
}

void main(){
    float w = texture(DiffMap3, vUV).r;

    vec3 albedo = mix(texture(DiffMap0, vUV).rgb,
                      texture(DiffMap2, vUV).rgb, w);
    float specR = mix(texture(SpecMap0, vUV).r,
                      texture(SpecMap1, vUV).r, w);

    vec3 N = normalize(vNormalWS);
    vec3 T = normalize(vTangentWS);
    vec3 B = normalize(vBinormalWS);

    if(twoSided == 1 && !gl_FrontFacing) { N = -N; B = -B; }

    mat3 TBN = mat3(T, B, N);
    vec3 tn0 = (isFlat0 == 1) ? vec3(0,0,1) : SampleNormalDXT5(texture(BumpMap0, vUV));
    vec3 tn1 = (isFlat1 == 1) ? vec3(0,0,1) : SampleNormalDXT5(texture(BumpMap1, vUV));
    vec3 tn  = BlendNormals(tn0, tn1, w);
    vec3 Nws = normalize(TBN * tn);

    gPosition   = vec4(vWorldPos, 1.0);
    gNormal     = vec4(Nws * 0.5 + 0.5, 1.0);
    gAlbedoSpec = vec4(albedo, specR);
}
)";


    const char* environment_fs = R"(
#version 330 core
uniform sampler2D DiffMap;
uniform sampler2D SpecMap;
uniform sampler2D BumpMap;
uniform sampler2D CubeMap0;

in vec3 vWorldPos;
in vec3 vNormalWS;
in vec2 vUV;
in vec3 vTangentWS;
in vec3 vBinormalWS;

layout(location=0) out vec4 gPosition;
layout(location=1) out vec4 gNormal;
layout(location=2) out vec4 gAlbedoSpec;

vec3 SampleNormalDXT5(vec4 s){ vec2 xy = s.wy*2.0-1.0; float z = sqrt(max(0.0,1.0-dot(xy,xy))); return vec3(xy,z); }
vec2 DirToLatLong(vec3 r){
    float u = atan(r.z, r.x) / (2.0*3.14159265359) + 0.5;
    float v = 0.5 - asin(clamp(r.y,-1.0,1.0)) / 3.14159265359;
    return vec2(u,v);
}

void main(){
    vec4 d = texture(DiffMap, vUV);
    vec4 s = texture(SpecMap, vUV);
    vec4 b = texture(BumpMap, vUV);

    vec3 N = normalize(vNormalWS);
    vec3 T = normalize(vTangentWS);
    vec3 B = normalize(vBinormalWS);
    mat3 TBN = mat3(T,B,N);
    vec3 Nws = normalize(TBN * SampleNormalDXT5(b));

    vec3 V = normalize(vWorldPos - vWorldPos);
    vec2 envUV = DirToLatLong(reflect(vec3(0,0,1), Nws));

    vec3 envCol = texture(CubeMap0, envUV).rgb;
    gPosition   = vec4(vWorldPos, 1.0);
    gNormal     = vec4(Nws*0.5+0.5, 1.0);
    gAlbedoSpec = vec4(mix(d.rgb, envCol, clamp(s.r*0.5,0.0,1.0)), s.r);
}

)";

    shaderSimpleLayer = make(deferred_vs, simplelayer_gbuffer_fs);
    shaderStandard = make(deferred_vs, deferred_fs);
    shaderLighting = make(lighting_vs, lighting_fs);
    shaderEnvironment = make(deferred_vs, environment_fs);
}


static inline glm::mat4 TRS_dx_to_gl(const std::array<float, 4>& P,
    const std::array<float, 4>& Q,
    const std::array<float, 4>& S)
{
    glm::vec3 p(P[0], P[1], -P[2]);
    glm::quat q(Q[3], -Q[0], -Q[1], Q[2]);
    if (glm::length(q) < 1e-6f) q = glm::quat(1, 0, 0, 0);
    q = glm::normalize(q);
    return glm::translate(glm::mat4(1.f), p) * glm::mat4_cast(q)
        * glm::scale(glm::mat4(1.f), glm::vec3(S[0], S[1], S[2]));
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

bool debug_output = false;

void InitDeferred(int w, int h) {
    glGenFramebuffers(1, &gBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);

    glGenTextures(1, &gPosition);
    glBindTexture(GL_TEXTURE_2D, gPosition);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gPosition, 0);

    glGenTextures(1, &gNormal);
    glBindTexture(GL_TEXTURE_2D, gNormal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gNormal, 0);

    glGenTextures(1, &gAlbedoSpec);
    glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gAlbedoSpec, 0);

    glGenTextures(1, &gDepth);
    glBindTexture(GL_TEXTURE_2D, gDepth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, w, h, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, gDepth, 0);

    GLenum attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
    glDrawBuffers(3, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) errout << "[GBuffer] fail\n";
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

inline void ResizeDeferred(int w, int h) {
    glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);

    if (gPosition) glDeleteTextures(1, &gPosition);
    if (gNormal) glDeleteTextures(1, &gNormal);
    if (gAlbedoSpec) glDeleteTextures(1, &gAlbedoSpec);
    if (gDepth) glDeleteTextures(1, &gDepth);

    glGenTextures(1, &gPosition);
    glBindTexture(GL_TEXTURE_2D, gPosition);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gPosition, 0);

    glGenTextures(1, &gNormal);
    glBindTexture(GL_TEXTURE_2D, gNormal);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gNormal, 0);

    glGenTextures(1, &gAlbedoSpec);
    glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, gAlbedoSpec, 0);

    glGenTextures(1, &gDepth);
    glBindTexture(GL_TEXTURE_2D, gDepth);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, w, h, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, gDepth, 0);

    GLenum attachments[3] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2 };
    glDrawBuffers(3, attachments);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        errout << "[GBuffer] Resize failed\n";

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}


GLuint LoadDDS(const std::string& path) {
    auto itc = gTexCache.find(path);
    if (itc != gTexCache.end()) return itc->second;

    std::string norm = path;
    if (norm.rfind("tex:", 0) == 0) norm = norm.substr(4);
    for (auto& ch : norm) if (ch == '/') ch = '\\';

    DirectX::ScratchImage image;
    DirectX::TexMetadata meta;
    HRESULT hr = DirectX::LoadFromDDSFile(
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

    if (debug_output) {
        std::cout << "[DDS INFO] " << norm
            << " | IsCubemap=" << (meta.IsCubemap() ? "true" : "false")
            << " | ArraySize=" << meta.arraySize
            << " | Mips=" << meta.mipLevels
            << " | Format=" << meta.format
            << " | " << meta.width << "x" << meta.height << "\n";
    }

    if (meta.IsCubemap()) {
        GLenum glFormat = 0;
        switch (meta.format) {
        case DXGI_FORMAT_BC1_UNORM: glFormat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT; break;
        case DXGI_FORMAT_BC2_UNORM: glFormat = GL_COMPRESSED_RGBA_S3TC_DXT3_EXT; break;
        case DXGI_FORMAT_BC3_UNORM: glFormat = GL_COMPRESSED_RGBA_S3TC_DXT5_EXT; break;
        case DXGI_FORMAT_BC1_UNORM_SRGB: glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT1_EXT; break;
        case DXGI_FORMAT_BC2_UNORM_SRGB: glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT3_EXT; break;
        case DXGI_FORMAT_BC3_UNORM_SRGB: glFormat = GL_COMPRESSED_SRGB_ALPHA_S3TC_DXT5_EXT; break;
        default: glFormat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT; break;
        }

        GLuint tex = 0;
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_CUBE_MAP, tex);

        for (size_t face = 0; face < 6; ++face) {
            for (size_t level = 0; level < meta.mipLevels; ++level) {
                const DirectX::Image* mip = image.GetImage(level, face, 0);
                glCompressedTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + (GLenum)face,
                    (GLint)level, glFormat,
                    (GLsizei)mip->width, (GLsizei)mip->height, 0,
                    (GLsizei)mip->slicePitch, mip->pixels);
            }
        }

        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

        glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
        gTexCache[path] = tex;
        return tex;
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
    default: glFormat = 0; break;
    }

    if (glFormat == 0) {
        DirectX::ScratchImage conv;
        const DirectX::Image* src = img;
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
            (GLsizei)src->width, (GLsizei)src->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, src->pixels);
        glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glBindTexture(GL_TEXTURE_2D, 0);

        gTexCache[path] = tex;
        return tex;
    }
    else {
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
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

        glBindTexture(GL_TEXTURE_2D, 0);
        gTexCache[path] = tex;
        return tex;
    }
}

void initFallbackTextures() {
    glGenTextures(1, &gWhiteTex);
    glBindTexture(GL_TEXTURE_2D, gWhiteTex);
    unsigned char white[4] = { 255, 255, 255, 255 };
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, white);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glGenTextures(1, &gBlackTex);
    glBindTexture(GL_TEXTURE_2D, gBlackTex);
    unsigned char black[4] = { 0, 0, 0, 255 };
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, black);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glGenTextures(1, &gFlatNormalTex);
    glBindTexture(GL_TEXTURE_2D, gFlatNormalTex);
    unsigned char flat[4] = { 128, 128, 255, 255 };
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, 1, 1, 0, GL_RGBA, GL_UNSIGNED_BYTE, flat);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glBindTexture(GL_TEXTURE_2D, 0);
}
struct Reporter {
    enum Type { Info, Warn, Err };
    std::string currentFile;

    void report(Type t, const std::string& msg) {
        const char* k = t == Info ? "INFO" : t == Warn ? "WARNING" : "ERROR";
        errout << "[" << k << "] ";
        if (!currentFile.empty()) {
            errout << "File '" << currentFile << "': ";
        }
        errout << msg << "\n";
    }
};

class Parser {
public:
    inline static const std::unordered_set<std::string> kAllowed = {
        "NEB3","3BEN",
        ">MDL","<MDL",">MND","<MND","LDM>","LDM<","DNM>","DNM<",
        "EOF_","_FOE","MODL","CHRN","CHSN","NRHC",
        "NFRT","TRFN","TNIS","XOBL","NSHC","LDOM",
        "MESH","HSEM","SMES","PGRI","IRGP","GRFS","STXT",
        "MATE","RDHS", "SEMF", "SPND", "HCRH",
        "POSI","ROTN","SCAL","LBOX","MNTP","SSTA",
        "TXTS","SFLT","TLFS","SFV4","CEVS","STUS","SCVU",
        "NSKF","TNJN","TNOJ",
        "SFRG","PEDU","PLPE","PACD","PROF","PBBO",
        "PRMN","PRMX","PGRV","PSTC","PTTX","PSTS","PVRM","PRRM",
        "PSRM","PPCT","PRRD","PSDL","PVAF","PDEL",
        "STUS","SSPI","SSTA",
        "HSAC","PTNM","SEMK","SEMJ","SEMU","IRGP","FKSN",
        "ISOP","NTOR","LACS","ETAM","SUTS","UVCS",
        "PSND",
        "BASE","SLPT","ANNO","SANI","ADDK","ADEK","ADSK",
        "QRFE","TFLP","NMSP","XMSP","LVSP","LVRP","EZSP","PLAP","SSMP",
        "NMTP","FLVP","RIAP","DERP","NRGP","ULBP",
        "XTTP","STSP","MRVP","MRRP","MRSP",
        "UDEP","TCPP","DCAP","VRGP","CTSP","FAVP","LEDP",
        "EPLP","CSLP","FORP","OBBP","DRRP","LDSP",
        "NMRP","XMRP",
        "EFRQ","PLFT","PSMN","PSMX","PSVL","PRVL","PSZE","PMSS",
        "PTMN","PVLF","PAIR","PRED","PGRN","PBLU","PALP",
        "MINA","TRAV","LNKS","LKSN","SKNL","NSKL"
    };


    Parser(Reporter& r, const Options& opt) :rep(r), options(opt) {}

    const std::vector<Node*>& getNodes() const { return n3node_list; }

    bool parse_file(const std::string& filepath) {
        this->filepath = filepath;
        rep.currentFile = filepath;
        std::ifstream f(filepath, std::ios::binary);
        if (!f) { rep.report(Reporter::Err, "Cannot open file"); return false; }
        if (!read_header(f)) return false;

        if (!(n3version == 1 || n3version == 2)) {
            if (options.ignore_version) rep.report(Reporter::Warn, "Unsupported version '" + std::to_string(n3version) + "'");
            else { rep.report(Reporter::Err, "Unsupported version '" + std::to_string(n3version) + "'"); return false; }
        }

        nodeStack.clear();
        currentNFRT.clear();
        std::vector<Node*> nodeStackObj;
        Node* animOwner = nullptr;

        bool done = false;
        while (!done && !f.eof()) {
            std::streampos tagPos = f.tellg();
            std::string tag; if (!readFourCC(f, tag)) break;

            if (tag == ">MDL" || tag == "LDM>") {
                std::string t; if (!readFourCC(f, t)) return false;
                n3modeltype = t;
                n3modelname = readString(f);
                std::cout << "LDM> (Begin Model Block)\n";
                std::cout << "  LDOM Type: " << t << "\n";
                std::cout << "  LDOM Model Name: " << n3modelname << "\n";
                continue;
            }
            else if (tag == "<MDL" || tag == "LDM<") {
                if (animOwner) {
                    std::cout << std::string((nodeDepth(animOwner) + 1) * 2, ' ') << "=== ANIMATION HIERARCHY END ===\n";
                    animOwner = nullptr;
                }
                std::cout << "LDM< (End Model Block)\n";
                done = true;
                continue;
            }
            else if (tag == "EOF_") {
                if (animOwner) {
                    std::cout << std::string((nodeDepth(animOwner) + 1) * 2, ' ') << "=== ANIMATION HIERARCHY END ===\n";
                    animOwner = nullptr;
                }
                std::cout << "_FOE (End of File)\n";
                done = true;
                continue;
            }
            else if (tag == ">MND" || tag == "DNM>") {
                std::string t; if (!readFourCC(f, t)) return false;
                std::string node_name = readString(f);

                auto new_node = std::make_unique<Node>();
                Node* n = new_node.get();
                if (!n) return false;

                n->node_type = t;
                n->node_name = node_name;
                n->node_parent = nodeStackObj.empty() ? nullptr : nodeStackObj.back();

                if (n->node_parent) n->node_parent->node_children.push_back(n);
                n3node_list.push_back(n);
                n3node_storage.push_back(std::move(new_node));
                nodeStackObj.push_back(n);

                std::cout << std::string(nodeStackObj.size() * 2, ' ') << "DNM> " << t << " \"" << node_name << "\"\n";

                if (t == "NRHC") {
                    logout << std::string((nodeStackObj.size() + 1) * 2, ' ') << "=== ANIMATION HIERARCHY START ===\n";
                    animOwner = n;
                }

                if (t == "NFRT") {
                    nodeStack.push_back(node_name);
                    currentNFRT = node_name;
                }

                if (t == "NSHC" || t == "SPND" || t == "SEMF") {
                    logout << std::string((nodeStackObj.size() + 1) * 2, ' ') << "[NODE TYPE] " << t << "\n";
                }

                continue;
            }
            else if (tag == "<MND" || tag == "DNM<") {
                if (!nodeStackObj.empty()) {
                    Node* closing = nodeStackObj.back();
                    nodeStackObj.pop_back();
                    if (closing) {
                        if (animOwner == closing) {
                            logout << std::string((nodeStackObj.size() + 2) * 2, ' ') << "=== ANIMATION HIERARCHY END ===\n";
                            animOwner = nullptr;
                        }

                        std::cout << std::string((nodeStackObj.size() + 1) * 2, ' ') << "DNM< (End " << closing->node_type << " \"" << closing->node_name << "\")\n";

                        if (closing->node_type == "NFRT") {
                            if (!nodeStack.empty()) nodeStack.pop_back();
                            currentNFRT = nodeStack.empty() ? std::string() : nodeStack.back();
                        }
                    }
                }
                continue;
            }
            else {
                if (!nodeStackObj.empty()) {
                    Node* current = nodeStackObj.back();
                    if (current && parse_node_tag(f, tag, *current)) continue;
                }
                std::streampos beforeSkip = f.tellg();
                f.seekg(tagPos);
                if (!skip_unknown(f)) continue;
                if (f.tellg() != beforeSkip) {
                    errout << "[SKIP] Unknown tag '" << tag << "' at offset " << tagPos << "\n";
                }
            }
        }

        dumpHierarchy();
        return true;
    }


    const std::vector<Node*>& nodes() const { return n3node_list; }
    int version() const { return n3version; }
    const std::string& model_type() const { return n3modeltype; }
    const std::string& model_name() const { return n3modelname; }

private:
    void printHierarchy(const Node* node, int depth = 0) const {
        for (int i = 0; i < depth; i++) std::cout << "    ";
        std::cout << node->node_type << " - " << node->node_name << "\n";
        for (auto* child : node->node_children) {
            printHierarchy(child, depth + 1);
        }
    }

    void dumpHierarchy() const {
        std::cout << "\n=== HIERARCHY TREE ===\n";
        for (auto* n : n3node_list) {
            if (!n->node_parent) {
                printHierarchy(n, 0);
            }
        }
        std::cout << "======================\n";

    }


    bool skip_unknown(std::ifstream& f) {
        std::streampos start = f.tellg();
        std::streampos endpos;
        { auto save = f.tellg(); f.seekg(0, std::ios::end); endpos = f.tellg(); f.seekg(save); }
        const int64_t MAX_SCAN = 4 * 1024 * 1024;
        char buf[4];
        for (int64_t i = 1; i < MAX_SCAN && (start + std::streamoff(i)) <= endpos - std::streamoff(4); ++i) {
            f.seekg(start + std::streamoff(i));
            if (!f.read(buf, 4)) return true;
            std::string t(buf, 4);
            if (is_valid_fourcc(t.c_str()) && kAllowed.find(t) != kAllowed.end()) {
                f.seekg(start + std::streamoff(i));
                return true;
            }
        }
        f.seekg(endpos);
        return true;
    }

    Reporter& rep;
    Options options;
    std::string filepath;
    int n3version = 0;
    std::string n3modeltype;
    std::string n3modelname;
    std::vector<std::unique_ptr<Node>> n3node_storage;
    std::vector<Node*> n3node_list;
    bool have_pending = false;
    std::vector<std::string> nodeStack;
    std::string currentNFRT;
    std::string lastNFRT;
    bool swapDataBytes = false;
    bool animBlockActive = false;
    int animDepth = -1;
    const Node* animOwner = nullptr;

    static uint16_t bswap16(uint16_t v) { return uint16_t((v >> 8) | (v << 8)); }
    static uint32_t bswap32(uint32_t v) { return (v >> 24) | ((v >> 8) & 0x0000FF00u) | ((v << 8) & 0x00FF0000u) | (v << 24); }

    bool readU16(std::ifstream& f, uint16_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 2)) return false; if (swapDataBytes) out = bswap16(out); return true; }
    bool readI32(std::ifstream& f, int32_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 4)) return false; if (swapDataBytes) out = (int32_t)bswap32((uint32_t)out); return true; }
    bool readU32(std::ifstream& f, uint32_t& out) { if (!f.read(reinterpret_cast<char*>(&out), 4)) return false; if (swapDataBytes) out = bswap32(out); return true; }
    bool readF32(std::ifstream& f, float& out) { uint32_t u; if (!readU32(f, u)) return false; std::memcpy(&out, &u, 4); return true; }

    bool readI8(std::ifstream& f, int8_t& out) {
        uint8_t c;
        if (!f.read(reinterpret_cast<char*>(&c), 1)) return false;
        out = static_cast<int8_t>(c);
        return true;
    }

    int nodeDepth(const Node* n) const {
        int d = 0;
        for (auto* p = n; p; p = p->node_parent) ++d;
        return d;
    }

    std::string indent(int d) const {
        if (d < 0) d = 0;
        return std::string(static_cast<size_t>(d) * 4, ' ');
    }


    bool readFourCC(std::ifstream& f, std::string& s) {
        char b[4];
        if (!f.read(b, 4)) return false;
        s.assign(b, 4);
        return true;
    }

    static inline bool is_valid_fourcc(const char* b) {
        for (int i = 0; i < 4; ++i) {
            unsigned char c = (unsigned char)b[i];
            if (!((c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_' || c == '<' || c == '>'))
                return false;
        }
        return true;
    }


    std::string readString(std::ifstream& f) {
        uint16_t n = 0; if (!readU16(f, n)) return {};
        std::string s; s.resize(n);
        if (n > 0) f.read(&s[0], n);
        return s;
    }

    bool read_header(std::ifstream& f) {
        char magic[4]; if (!f.read(magic, 4)) { rep.report(Reporter::Err, "Failed reading header"); return false; }
        if (std::memcmp(magic, "NEB3", 4) != 0 && std::memcmp(magic, "3BEN", 4) != 0) {
            rep.report(Reporter::Err, std::string("Invalid file, unknown fourCC '") + std::string(magic, 4) + "'");
            return false;
        }

        uint32_t rawVer = 0;
        if (!f.read(reinterpret_cast<char*>(&rawVer), 4)) return false;
        uint32_t le = rawVer, be = bswap32(rawVer);
        if (le == 1 || le == 2) { swapDataBytes = false; n3version = (int)le; }
        else if (be == 1 || be == 2) { swapDataBytes = true; n3version = (int)be; }
        else { rep.report(Reporter::Err, "Unsupported version field"); return false; }

        return true;
    }

    bool parse_node_tag(std::ifstream& f, const std::string& tag, Node& node) {
        // NSHC, SEMF are node types (handled by DNM>), not tags

        if (tag == ">MND" || tag == "DNM>") return true;
        if (tag == "<MND" || tag == "DNM<") return true;

        if (tag == "HSEM") return parse_mesh(f, node);
        if (tag == "IRGP") return parse_pgri(f, node);
        if (tag == "POSI") return parse_posi(f, node);
        if (tag == "ROTN") return parse_rotn(f, node);
        if (tag == "SCAL") return parse_scal(f, node);
        if (tag == "TXTS") return parse_stxt(f, node);
        if (tag == "RDHS") return parse_shdr(f, node);
        if (tag == "CEVS") return parse_cevs(f, node);
        if (tag == "XOBL") return parse_xobl(f, node);
        if (tag == "TLFS") return parse_tlfs(f);
        if (tag == "TNIS") return parse_tnis(f);
        if (tag == "PTNM") return parse_mntp(f, node);
        if (tag == "HSAC") return parse_hsac(f, node);
        if (tag == "NRHC") return parse_nrhc(f, node);
        if (tag == "MINA") return parse_mina(f, node);
        if (tag == "TRAV") return parse_trav(f, node);
        if (tag == "TNJN") return parse_tnjn(f, node);
        if (tag == "LKSN") return parse_lksn(f, node);
        if (tag == "LNKS") return parse_lnks(f, node);
        if (tag == "HCRH") return parse_hcrh(f, node);

        return false;
    }
    bool parse_hcrh(std::ifstream& f, Node& node) {
        int8_t val; if (!readI8(f, val)) return false;
        std::cout << std::string((nodeDepth(&node) + 1) * 2, ' ') << "HCRH val=" << (int)val << "\n";
        return true;
    }
    bool parse_hsac(std::ifstream& f, Node& node) {
        int8_t mode; if (!readI8(f, mode)) return false;
        std::cout << std::string((nodeDepth(&node) + 1) * 2, ' ') << "HSAC" << " mode=" << (int)mode << "\n";
        return true;
    }

    bool parse_nrhc(std::ifstream& f, Node& node) {
        node.node_type = "CHRN";
        // no need to read, we already do that in parse_file
        return true;
    }
    bool parse_mina(std::ifstream& f, Node& node)
    {
        std::string animRes = readString(f);
        std::cout << std::string((nodeDepth(&node) + 1) * 2, ' ') << "MINA" << " \"" << animRes << "\"\n";
        return true;
    }

    bool parse_trav(std::ifstream& f, Node& node)
    {
        std::string varRes = readString(f);
        std::cout << std::string((nodeDepth(&node) + 1) * 2, ' ') << "TRAV" << " \"" << varRes << "\"\n";
        return true;
    }
    bool parse_tnjn(std::ifstream& f, Node& node)
    {
        int32_t numJoints; if (!readI32(f, numJoints)) return false;
        std::cout << std::string((nodeDepth(&node) + 1) * 2, ' ') << "TNJN" << " count=" << numJoints << "\n";
        return true;
    }

    bool parse_lksn(std::ifstream& f, Node& node)
    {
        int32_t count;
        if (!readI32(f, count)) return false;
        std::cout << "[LKSN] NumSkinLists = " << count << "\n";
        return true;
    }

    bool parse_lnks(std::ifstream& f, Node& node) {
        std::string skinListName = readString(f);
        int32_t num = 0;
        if (!readI32(f, num)) return false;
        std::cout << "[LNKS] SkinList '" << skinListName << "' count=" << num << "\n";

        for (int i = 0; i < num; i++) {
            std::string skin = readString(f);
            std::cout << "   skin[" << i << "] = " << skin << "\n";
        }

        int8_t unknownByte;
        if (!readI8(f, unknownByte)) return false;

        std::string varName = readString(f);
        if (!varName.empty()) {
            std::cout << "   [VARIATION] " << varName << "\n";
        }
        return true;
    }

    bool parse_mntp(std::ifstream& f, Node& node) {
        std::string s = readString(f);
        node.mat = s;
        std::cout << "        MNTP " << s << "\n";
        return true;
    }

    bool parse_cevs(std::ifstream& f, Node& node) {
        uint16_t nameLen;
        if (!readU16(f, nameLen)) return false;
        std::string param(nameLen, '\0');
        if (!f.read(param.data(), nameLen)) return false;

        std::array<float, 4> v{};
        for (int i = 0; i < 4; i++)
            if (!readF32(f, v[i])) return false;

        std::cout << "[CEVS] " << node.node_name << " param=" << param
            << " val=(" << v[0] << "," << v[1] << "," << v[2] << "," << v[3] << ")\n";
        return true;
    }

    bool parse_xobl(std::ifstream& f, Node& node) {
        float data[8];
        for (int i = 0; i < 8; i++) {
            if (!readF32(f, data[i])) return false;
        }
        std::cout << "[XOBL] node=" << node.node_name
            << " pos(" << data[0] << "," << data[1] << "," << data[2] << ")"
            << " rot(" << data[3] << "," << data[4] << "," << data[5] << "," << data[6] << ")"
            << " scl(" << data[7] << "," << data[7] << "," << data[7] << ")\n";
        return true;
    }

    bool parse_shdr(std::ifstream& f, Node& node) {
        node.shdr = readString(f);
        while (!node.shdr.empty() && (node.shdr.back() <= 32 || node.shdr.back() == '\0'))
            node.shdr.pop_back();
        logout << "        SHDR = " << node.shdr << "\n";
        return true;
    }

    bool parse_tlfs(std::ifstream& f) {
        std::string key = readString(f);
        float val = 0.f; if (!readF32(f, val)) return false;
        std::cout << "        TLFS " << key << " = " << val << "\n";
        return true;
    }

    bool parse_tnis(std::ifstream& f) {
        std::string key = readString(f);
        int32_t v = 0; if (!readI32(f, v)) return false;
        std::cout << "        TNIS " << key << " = " << v << "\n";
        return true;
    }

    static inline bool n_isdenormal(float s) {
        uint32_t u = *(const uint32_t*)&s;
        return (u & 0x7F800000u) == 0;
    }

    static inline float n_undenormalize(float f) {
        if (n_isdenormal(f) || std::fabs(f) < 1.0e-30f)
            return 0.0f;
        return f;
    }

    bool parse_mesh(std::ifstream& f, Node& node) {
        node.mesh_ressource_id = readString(f);
        logout << "        MESH " << node.mesh_ressource_id << "\n";
        return true;
    }
    bool parse_pgri(std::ifstream& f, Node& node) {
        int32_t v; if (!readI32(f, v)) return false;
        node.primitive_group_idx = v;
        std::cout << "        PGRI " << node.primitive_group_idx << "\n";
        return true;
    }

    bool parse_posi(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++)
            if (!readF32(f, node.position[i])) return false;

        logout << "[POSI] "
            << "x=" << node.position[0] << " "
            << "y=" << node.position[1] << " "
            << "z=" << node.position[2] << " "
            << "w=" << node.position[3] << "\n";

        return true;
    }

    bool parse_rotn(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++)
            if (!readF32(f, node.rotation[i])) return false;

        glm::vec4 v(node.rotation[0], node.rotation[1], node.rotation[2], node.rotation[3]);
        std::cout << "ROTN " << glm::to_string(v) << "\n";
        return true;
    }

    bool parse_scal(std::ifstream& f, Node& node) {
        for (int i = 0; i < 4; i++)
            if (!readF32(f, node.scale[i])) return false;

        glm::vec4 v(node.scale[0], node.scale[1], node.scale[2], node.scale[3]);
        std::cout << "SCAL " << glm::to_string(v) << "\n";
        return true;
    }

    bool parse_stxt(std::ifstream& f, Node& node) {
        std::string tex_type = readString(f);
        std::string tex_name = readString(f);
        node.shader_textures[tex_type] = tex_name;
        std::cout << "        STXT " << tex_type << " => " << tex_name << "\n";
        return true;
    }
};

std::vector<DrawCmd> BuildDrawsWithTransform(const Parser& parser, const glm::vec3& pos, const glm::quat& rot, const glm::vec3& scale);

GLuint LoadNebTexSmart(const std::string& path) {
    if (path.find("system/white.dds") != std::string::npos) return gWhiteTex;
    if (path.find("system/black.dds") != std::string::npos) return gBlackTex;
    if (path.find("system/nobump.dds") != std::string::npos) return gFlatNormalTex;

    GLuint tex = LoadDDS(path);
    return tex;
}

static GLuint loadNebTex(const std::unordered_map<std::string, std::string>& m, const char* key) {
    auto it = m.find(key);
    if (it == m.end()) return 0;
    std::string t = it->second;
    if (t.rfind("tex:", 0) == 0) t = t.substr(4);

    return LoadNebTexSmart(std::string("C:/drasa_online/work/textures/") + t + ".dds");
}

void InitSamplers() {
    glGenSamplers(1, &gSamplerRepeat);
    glSamplerParameteri(gSamplerRepeat, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glSamplerParameteri(gSamplerRepeat, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glSamplerParameteri(gSamplerRepeat, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glSamplerParameteri(gSamplerRepeat, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glGenSamplers(1, &gSamplerClamp);
    glSamplerParameteri(gSamplerClamp, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glSamplerParameteri(gSamplerClamp, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glSamplerParameteri(gSamplerClamp, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glSamplerParameteri(gSamplerClamp, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
}


void appendN3WTransform(const std::string& path,
    const glm::vec3& pos, const glm::quat& rot, const glm::vec3& scale)
{
    Reporter rep;
    Options opt;
    auto p = std::make_unique<Parser>(rep, opt);
    if (!p->parse_file(path)) {
        errout << "[FAIL] parse " << path << "\n";
        return;
    }

    auto draws = BuildDrawsWithTransform(*p, pos, rot, scale);
    for (auto& dc : draws)
        withTransform.push_back(std::move(dc));
}

void InitScreenQuad() {
    float quadVertices[] = {
        -1.0f,-1.0f, 0.0f,0.0f,
         1.0f,-1.0f, 1.0f,0.0f,
        -1.0f, 1.0f, 0.0f,1.0f,
         1.0f, 1.0f, 1.0f,1.0f
    };
    glGenVertexArrays(1, &quadVAO);
    glGenBuffers(1, &quadVBO);
    glBindVertexArray(quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);
}
void RenderLightingPass(GLuint shaderLighting, GLuint quadVAO, const glm::vec3& camPos) {
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glDisable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT);
    glUseProgram(shaderLighting);

    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, gPosition);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, gNormal);
    glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, gAlbedoSpec);
    glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, gDepth);

    glUniform1i(glGetUniformLocation(shaderLighting, "gPosition"), 0);
    glUniform1i(glGetUniformLocation(shaderLighting, "gNormal"), 1);
    glUniform1i(glGetUniformLocation(shaderLighting, "gAlbedoSpec"), 2);
    glUniform1i(glGetUniformLocation(shaderLighting, "gDepth"), 3);
    glUniform3fv(glGetUniformLocation(shaderLighting, "camPos"), 1, glm::value_ptr(camPos));

    GLint linked = 0;
    glGetProgramiv(shaderLighting, GL_LINK_STATUS, &linked);
    if (!linked) {
        std::cerr << "[LIGHT] shaderLighting link failed\n";
        return;
    }

    GLint vaoValid = 0;
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &vaoValid);
    if (!vaoValid) {
        glBindVertexArray(quadVAO);
    }

    GLint boundTex = 0;
    glActiveTexture(GL_TEXTURE2);
    glGetIntegerv(GL_TEXTURE_BINDING_2D, &boundTex);

    if (boundTex == 0) {
        std::cerr << "[LIGHT] gAlbedoSpec not bound\n";
        return;
    }

    int w = 0, h = 0;
    glBindTexture(GL_TEXTURE_2D, boundTex);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);

    if (w == 0 || h == 0) {
        std::cerr << "[LIGHT] invalid gAlbedoSpec dimensions\n";
        return;
    }

    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}

int main(int argc, char* argv[])
{
    try {
        InitGLFW();
        InitSHDR();
        initFallbackTextures();
        InitSamplers();

        double lastFrame = 0.0;
        double deltaTime = 0.0;

        glfwSetWindowUserPointer(window, &camera);

        Reporter rep; Options opt; Parser parser(rep, opt);

        // e_z0110_dun viscanium
        // e_z0120_dun some arena
        // e_b1205_dun_chr_2012_01 xmas map
        // e_nye_dungeon_01 January event
        // f4103_atf old 6v6 pvp map

        appendN3WTransform("C:\\drasa_online\\work\\models\\t001_hub\\deco_box_02.n3",
            { 0,0,0 }, { 0,0,0,0 }, { 1,1,1 });

        std::stable_sort(withTransform.begin(), withTransform.end(),
            [](const DrawCmd& a, const DrawCmd& b) { return a.shdr < b.shdr; });

        glfwSetWindowUserPointer(window, &camera);
        glfwSetCursorPosCallback(window, Camera::mouse_callback);
        glfwSetScrollCallback(window, Camera::scroll_callback);
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        InitScreenQuad();

        glm::vec3 camPos = camera.getPosition();
        // std::cout << "[CAMERA] pos=(" << camPos.x << "," << camPos.y << "," << camPos.z << ")\n";

        while (!glfwWindowShouldClose(window))
        {
            const double currentFrame = glfwGetTime();
            deltaTime = currentFrame - lastFrame;
            lastFrame = currentFrame;

            camera.handleInput(window, static_cast<float>(deltaTime));
            glm::vec3 cameraWorldPos = camera.getPosition();

            glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);
            glDepthFunc(GL_LESS);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            glDisable(GL_BLEND);

            for (auto& obj : withTransform)
            {
                if (obj.mesh.vao == 0 || obj.mesh.idx.empty()) {
                    std::cerr << "[RENDER] Skipping invalid mesh for " << obj.name << "\n";
                    continue;
                }

                GLuint prog = 0;
                if (obj.shdr == "shd:simplelayer")       prog = shaderStandard;
                else if (obj.shdr == "shd:environment")  prog = shaderEnvironment;
                else if (obj.shdr == "shd:decal")        prog = shaderStandard;
                else                                     prog = shaderStandard;

                glUseProgram(prog);

                glUniformMatrix4fv(glGetUniformLocation(prog, "projection"), 1, GL_FALSE,
                    glm::value_ptr(camera.getProjectionMatrix()));
                glUniformMatrix4fv(glGetUniformLocation(prog, "view"), 1, GL_FALSE,
                    glm::value_ptr(camera.getViewMatrix()));
                glUniformMatrix4fv(glGetUniformLocation(prog, "model"), 1, GL_FALSE,
                    glm::value_ptr(obj.worldMatrix));

                glUniform1i(glGetUniformLocation(prog, "alphaTest"),
                    (obj.mat == "AlphaTest") ? 1 : 0);
                glUniform1i(glGetUniformLocation(prog, "twoSided"),
                    (obj.mat == "AlphaTest" || obj.mat == "Decal") ? 1 : 0);
                glUniform1i(glGetUniformLocation(prog, "isFlatNormal"),
                    (obj.tex[2] == gFlatNormalTex) ? 1 : 0);
                glUniform1f(glGetUniformLocation(prog, "alphaCutoff"), 0.5f);

                if (obj.shdr == "shd:simplelayer")
                {
                    glUniform1i(glGetUniformLocation(prog, "DiffMap0"), 0);
                    glUniform1i(glGetUniformLocation(prog, "SpecMap0"), 1);
                    glUniform1i(glGetUniformLocation(prog, "BumpMap0"), 2);
                    glUniform1i(glGetUniformLocation(prog, "DiffMap2"), 4);
                    glUniform1i(glGetUniformLocation(prog, "SpecMap1"), 5);
                    glUniform1i(glGetUniformLocation(prog, "BumpMap1"), 6);
                    glUniform1i(glGetUniformLocation(prog, "DiffMap3"), 7);
                    glUniform1i(glGetUniformLocation(prog, "isFlat0"),
                        (obj.tex[2] == gFlatNormalTex) ? 1 : 0);
                    glUniform1i(glGetUniformLocation(prog, "isFlat1"),
                        (obj.tex[6] == gFlatNormalTex) ? 1 : 0);

                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[0] ? obj.tex[0] : gWhiteTex);
                    glBindSampler(0, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE1);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[1] ? obj.tex[1] : gWhiteTex);
                    glBindSampler(1, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE2);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[2] ? obj.tex[2] : gFlatNormalTex);
                    glBindSampler(2, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE4);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[4] ? obj.tex[4] : gWhiteTex);
                    glBindSampler(4, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE5);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[5] ? obj.tex[5] : gWhiteTex);
                    glBindSampler(5, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE6);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[6] ? obj.tex[6] : gFlatNormalTex);
                    glBindSampler(6, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE7);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[7] ? obj.tex[7] : gWhiteTex);
                    glBindSampler(7, gSamplerClamp);
                }
                else if (obj.shdr == "shd:environment")
                {
                    glUniform1i(glGetUniformLocation(prog, "DiffMap"), 0);
                    glUniform1i(glGetUniformLocation(prog, "SpecMap"), 1);
                    glUniform1i(glGetUniformLocation(prog, "BumpMap"), 2);
                    glUniform1i(glGetUniformLocation(prog, "CubeMap0"), 4);

                    glm::vec3 cp = camera.getPosition();
                    glUniform3fv(glGetUniformLocation(prog, "cameraPos"), 1, glm::value_ptr(cp));

                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[0] ? obj.tex[0] : gWhiteTex);
                    glBindSampler(0, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE1);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[1] ? obj.tex[1] : gWhiteTex);
                    glBindSampler(1, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE2);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[2] ? obj.tex[2] : gFlatNormalTex);
                    glBindSampler(2, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE4);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[4] ? obj.tex[4] : 0);
                    glBindSampler(4, gSamplerClamp);
                }
                else // shd:standard
                {
                    glUniform1i(glGetUniformLocation(prog, "DiffMap"), 0);
                    glUniform1i(glGetUniformLocation(prog, "SpecMap"), 1);
                    glUniform1i(glGetUniformLocation(prog, "BumpMap"), 2);
                    glUniform1i(glGetUniformLocation(prog, "EmsvMap"), 3);

                    glActiveTexture(GL_TEXTURE0);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[0] ? obj.tex[0] : gWhiteTex);
                    glBindSampler(0, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE1);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[1] ? obj.tex[1] : gWhiteTex);
                    glBindSampler(1, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE2);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[2] ? obj.tex[2] : gFlatNormalTex);
                    glBindSampler(2, gSamplerRepeat);

                    glActiveTexture(GL_TEXTURE3);
                    glBindTexture(GL_TEXTURE_2D, obj.tex[3] ? obj.tex[3] : gBlackTex);
                    glBindSampler(3, gSamplerRepeat);
                }

                glBindVertexArray(obj.mesh.vao);

                if (obj.group >= 0 && obj.group < (int)obj.mesh.groups.size()) {
                    const Nvx2Group& g = obj.mesh.groups[obj.group];
                    glDrawElements(GL_TRIANGLES, g.indexCount(), GL_UNSIGNED_INT,
                        (void*)(intptr_t)(g.firstIndex() * sizeof(uint32_t)));
                }
                else {
                    for (size_t gi = 0; gi < obj.mesh.groups.size(); ++gi) {
                        const Nvx2Group& g = obj.mesh.groups[gi];
                        glDrawElements(GL_TRIANGLES, g.indexCount(), GL_UNSIGNED_INT,
                            (void*)(intptr_t)(g.firstIndex() * sizeof(uint32_t)));
                    }
                }
            }

            RenderLightingPass(shaderLighting, quadVAO, cameraWorldPos);

            glDepthMask(GL_TRUE);
            glDepthFunc(GL_LESS);
            glEnable(GL_CULL_FACE);

            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        for (auto& obj : withTransform) {
            glDeleteVertexArrays(1, &obj.mesh.vao);
            glDeleteBuffers(1, &obj.mesh.vbo);
            glDeleteBuffers(1, &obj.mesh.ebo);
        }
        withTransform.clear();
        glDeleteProgram(shaderProgram);
        glfwTerminate();
        return 0;
    }
    catch (const std::exception& e) {
        errout << "[FATAL] " << e.what() << "\n";
        return 1;
    }
}

std::vector<DrawCmd> BuildDrawsWithTransform(const Parser& parser, const glm::vec3& pos, const glm::quat& rot, const glm::vec3& scale)
{
    size_t countBuilt = 0, countDecal = 0, countNoMesh = 0, countNoDiff = 0;
    std::vector<DrawCmd> out;

    glm::mat4 S = glm::scale(glm::mat4(1.0f), scale);
    glm::mat4 R = glm::mat4_cast(rot);
    glm::mat4 T = glm::translate(glm::mat4(1.0f), pos);
    glm::mat4 trs = T * R * S;

    glm::vec3 trsPos = glm::vec3(trs[3]);
    // std::cout << "[BUILD TRANSFORM] TRS matrix position=(" << trsPos.x << "," << trsPos.y << "," << trsPos.z << ")\n";

    std::string modelName = parser.model_name();
    std::transform(modelName.begin(), modelName.end(), modelName.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    for (auto* node : parser.getNodes()) {
        if (!node) continue;

        /*
        std::cout << "[NODE CHECK] type=" << node->node_type
            << " name='" << node->node_name << "'"
            << " mesh='" << node->mesh_ressource_id << "'"
            << " shdr='" << node->shdr << "'"
            << " depth=" << GetNodeDepth(node) << "\n";*/

        if (node->mesh_ressource_id.empty()) { countNoMesh++; continue; }

        /*
        std::cout << "[NODE CHECK] type=" << node->node_type
            << " name='" << node->node_name << "'"
            << " mesh='" << node->mesh_ressource_id << "'"
            << " shdr='" << node->shdr << "'"
            << " depth=" << GetNodeDepth(node) << "\n";*/

        // NORMAL NON-DECAL PROCESSING
        if (node->mesh_ressource_id.empty()) {
            std::cout << "  -> SKIPPED: no mesh_ressource_id\n";
            continue;
        }


        auto it = node->shader_textures.find("DiffMap0");
        if (it == node->shader_textures.end())
            it = node->shader_textures.find("DiffMap");

        // std::cout << "  [TEX CHECK] has DiffMap: " << (it != node->shader_textures.end()) << "\n";

        if (it == node->shader_textures.end()) {
            std::cout << "  -> SKIPPED: no DiffMap0/DiffMap\n";
            continue;
        }

        const std::string& d = it->second;
        std::cout << "  [TEX CHECK] DiffMap value: '" << d << "' len=" << d.size() << "\n";

        if (d.empty() || d.size() < 8) {
            std::cout << "  -> SKIPPED: DiffMap too short\n";
           continue;
        }
        if (d.find("system/white") != std::string::npos) {
            std::cout << "  -> SKIPPED: DiffMap is system/white\n";
           continue;
        }
        if (d.find("system/black") != std::string::npos) {
            std::cout << "  -> SKIPPED: DiffMap is system/black\n";
           continue;
        }
        if (d.find("system/nobump") != std::string::npos) {
            std::cout << "  -> SKIPPED: DiffMap is system/nobump\n";
           continue;
        }


        if (it == node->shader_textures.end()) countNoDiff++;

        // if (node->shdr.empty() && node->shader_textures.empty()) continue;

        bool hasAny = node->shader_textures.count("DiffMap0") || node->shader_textures.count("DiffMap")
            || node->shader_textures.count("SpecMap0") || node->shader_textures.count("SpecMap")
            || node->shader_textures.count("BumpMap0") || node->shader_textures.count("BumpMap")
            || node->shader_textures.count("EmsvMap0");
        if (!hasAny) continue;

        DrawCmd dc;
        dc.group = node->primitive_group_idx;
        dc.nodeName = node->node_name;
        dc.mat = node->mat;
        dc.shdr = node->shdr;
        dc.shader_textures = node->shader_textures;

        // logout << dc.shdr << " ";

        /*
        if (node->shdr == "shd:particle" ||
            node->shdr == "shd:uvanimated" ||
            node->shdr == "shd:refraction") {
            continue;
        }*/

        glm::mat4 local = TRS_dx_to_gl(node->position, node->rotation, node->scale);
        for (const Node* p = node->node_parent; p; p = p->node_parent)
            local = TRS_dx_to_gl(p->position, p->rotation, p->scale) * local;
        dc.worldMatrix = trs * local;

        std::string meshPath = node->mesh_ressource_id;
        if (meshPath.rfind("msh:", 0) == 0) meshPath = meshPath.substr(4);
        meshPath = MESHES_ROOT + meshPath;

        // std::cout << "  [MESH] Loading: " << meshPath << "\n";

        if (!LoadNVX2(meshPath, dc.mesh)) {
            errout << "  -> FAILED to load mesh!\n";
            continue;
        }

        /*
        std::cout << "  [MESH] Loaded successfully: verts=" << dc.mesh.verts.size()
            << " indices=" << dc.mesh.idx.size()
            << " groups=" << dc.mesh.groups.size() << "\n";*/

        auto pick = [&](const char* a, const char* b = "") {
            GLuint t = loadNebTex(node->shader_textures, a);
            if (t) return t;
            if (b[0]) return loadNebTex(node->shader_textures, b);
            return (GLuint)0;
            };

        if (node->shdr == "shd:simplelayer") {
            dc.tex[0] = pick("DiffMap0", "DiffMap");
            dc.tex[1] = pick("SpecMap0", "SpecMap");
            dc.tex[2] = pick("BumpMap0", "BumpMap");
            dc.tex[3] = pick("EmsvMap0", "EmsvMap");
            dc.tex[4] = pick("DiffMap2");
            dc.tex[5] = pick("SpecMap1");
            dc.tex[6] = pick("BumpMap1");
            dc.tex[7] = pick("DiffMap3");
            for (int i = 0; i < 8; ++i) {
                dc.has[i] = dc.tex[i] != 0;
                if (!dc.has[i]) dc.tex[i] = (i == 0 || i == 4) ? gWhiteTex : gBlackTex;
            }
        }
        else if (node->shdr == "shd:environment") {
            dc.tex[0] = pick("DiffMap0", "DiffMap");
            dc.tex[1] = pick("SpecMap0", "SpecMap");
            dc.tex[2] = pick("BumpMap0", "BumpMap");
            dc.tex[3] = pick("EmsvMap0", "EmsvMap");
            auto it = node->shader_textures.find("CubeMap0");
            if (it != node->shader_textures.end()) {
                std::string base = it->second;
                if (base.rfind("tex:", 0) == 0) base = base.substr(4);
                dc.tex[4] = LoadDDS(std::string("C:/drasa_online/work/textures/") + base);
            }
            else {
                dc.tex[4] = 0;
            }
            for (int i = 5; i < 8; ++i) dc.tex[i] = 0;
            for (int i = 0; i < 8; ++i) {
                dc.has[i] = (dc.tex[i] != 0 && dc.tex[i] != gBlackTex && dc.tex[i] != gWhiteTex);
                if (dc.tex[i] == 0) dc.tex[i] = (i == 0) ? gWhiteTex : gBlackTex;
            }
        }
        else {
            dc.tex[0] = pick("DiffMap0", "DiffMap");
            dc.tex[1] = pick("SpecMap0", "SpecMap");
            dc.tex[2] = pick("BumpMap0", "BumpMap");
            dc.tex[3] = pick("EmsvMap0", "EmsvMap");
            for (int i = 4; i < 8; ++i) dc.tex[i] = 0;
            for (int i = 0; i < 8; ++i) {
                dc.has[i] = (dc.tex[i] != 0 && dc.tex[i] != gBlackTex && dc.tex[i] != gWhiteTex);
                if (dc.tex[i] == 0) dc.tex[i] = (i == 0) ? gWhiteTex : gBlackTex;
            }
        }

        dc.name = node->mesh_ressource_id;
        out.push_back(std::move(dc)); countBuilt++;
    }
    return out;
}