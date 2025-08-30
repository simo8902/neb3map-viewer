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

#define NOMINMAX
#include <DirectXTex/DirectXTex.h>
#include <random>

namespace fs = std::filesystem;

const std::vector<int> SUPPORTED_VERSIONS = { 1, 2 };

static void CheckShader(GLuint sh, const char* name) {
    GLint ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        GLint len = 0; glGetShaderiv(sh, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0'); glGetShaderInfoLog(sh, len, nullptr, log.data());
        std::cerr << "[SHADER COMPILE ERROR] " << name << "\n" << log << "\n";
    }
}
static void CheckProgram(GLuint prog) {
    GLint ok = 0; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        GLint len = 0; glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0'); glGetProgramInfoLog(prog, len, nullptr, log.data());
        std::cerr << "[PROGRAM LINK ERROR]\n" << log << "\n";
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

bool debug_output = true;

GLuint LoadTexture2D(const std::string& path) {
    std::wstring wpath(path.begin(), path.end());

    DirectX::ScratchImage image;
    DirectX::TexMetadata meta;
    HRESULT hr = DirectX::LoadFromDDSFile(wpath.c_str(),
        DirectX::DDS_FLAGS_NONE,
        &meta,
        image);
    if (FAILED(hr)) {
        std::cerr << "[ERROR] Failed to load DDS: " << path << "\n";
        return 0;
    }

    const DirectX::Image* img = image.GetImage(0, 0, 0);
    if (!img) {
        std::cerr << "[ERROR] DDS file has no image: " << path << "\n";
        return 0;
    }

    GLenum glFormat = 0;
    switch (meta.format) {
    case DXGI_FORMAT_BC1_UNORM:
    case DXGI_FORMAT_BC1_UNORM_SRGB:
        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
        break;
    case DXGI_FORMAT_BC2_UNORM:
    case DXGI_FORMAT_BC2_UNORM_SRGB:
        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT3_EXT;
        break;
    case DXGI_FORMAT_BC3_UNORM:
    case DXGI_FORMAT_BC3_UNORM_SRGB:
        glFormat = GL_COMPRESSED_RGBA_S3TC_DXT5_EXT;
        break;
    default:
        // Fallback to RGBA8 if unsupported
        DirectX::ScratchImage conv;
        hr = DirectX::Decompress(*img, DXGI_FORMAT_R8G8B8A8_UNORM, conv);
        if (FAILED(hr)) {
            std::cerr << "[ERROR] Unsupported DDS format " << meta.format << " for " << path << "\n";
            return 0;
        }
        img = conv.GetImage(0, 0, 0);
        GLuint tex;
        glGenTextures(1, &tex);
        glBindTexture(GL_TEXTURE_2D, tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8,
            (GLsizei)img->width, (GLsizei)img->height,
            0, GL_RGBA, GL_UNSIGNED_BYTE, img->pixels);
        glGenerateMipmap(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, 0);
        return tex;
    }

    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);

    size_t offset = 0;
    for (size_t level = 0; level < meta.mipLevels; ++level) {
        const DirectX::Image* mip = image.GetImage(level, 0, 0);
        GLsizei w = static_cast<GLsizei>(mip->width);
        GLsizei h = static_cast<GLsizei>(mip->height);
        GLsizei size = static_cast<GLsizei>(mip->slicePitch);

        glCompressedTexImage2D(GL_TEXTURE_2D, (GLint)level,
            glFormat,
            w, h, 0,
            size, mip->pixels);
        offset += size;
    }

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
        meta.mipLevels > 1 ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glBindTexture(GL_TEXTURE_2D, 0);
    return tex;
}


struct Node {
    std::string node_name;
    std::string node_type;
    std::shared_ptr<Node> node_parent;
    std::vector<std::shared_ptr<Node>> node_children;

    std::unordered_map<std::string, std::string> shader_textures;
    std::unordered_map<std::string, float> shader_parameters;

    std::string mesh_resource_id;
    std::string model_node_type;
};

class Parser {
private:
    std::string filepath;
    std::ifstream n3file;
    std::string byteorder;
    int n3version = 0;
    std::string n3modeltype;
    std::string n3modelname;
    std::vector<std::shared_ptr<Node>> n3node_list;

public:
    Parser() {}

    const std::vector<std::shared_ptr<Node>>& getNodes() const { return n3node_list; }

    void report(const std::string& rep_type, const std::string& rep_msg) {
        std::cout << "[" << rep_type << "] " << rep_msg << std::endl;
    }

    template<typename T>
    T read_n3_value() {
        T value{};
        n3file.read(reinterpret_cast<char*>(&value), sizeof(T));
        return value;
    }

    std::string read_n3_string() {
        uint16_t str_len = read_n3_value<uint16_t>();
        std::vector<char> buffer(str_len);
        n3file.read(buffer.data(), str_len);
        return std::string(buffer.begin(), buffer.end());
    }

    std::string read_n3_fourcc() {
        uint32_t val;
        n3file.read(reinterpret_cast<char*>(&val), sizeof(val));

        char buf[4];
        if (byteorder == "little") {
            buf[0] = val & 0xFF;
            buf[1] = (val >> 8) & 0xFF;
            buf[2] = (val >> 16) & 0xFF;
            buf[3] = (val >> 24) & 0xFF;
        }
        else {
            buf[3] = val & 0xFF;
            buf[2] = (val >> 8) & 0xFF;
            buf[1] = (val >> 16) & 0xFF;
            buf[0] = (val >> 24) & 0xFF;
        }

        return std::string(buf, 4);
    }

    bool parse_tag_shape(const std::string& tag_4cc, std::shared_ptr<Node>& node) {
        if (tag_4cc == "MESH") {
            node->mesh_resource_id = read_n3_string();
            std::cout << "        mesh_res_id=" << node->mesh_resource_id << std::endl;
        }
        else {
            return false;
        }
        return true;
    }

    bool skip4cc(const std::string& tag_4cc, std::shared_ptr<Node>& node) {
        auto isAllowed = [](unsigned char c)->bool { return (c >= 'A' && c <= 'Z') || c == '<' || c == '>' || c == '_' || (c >= '0' && c <= '9'); };
        auto peekIsFourCC = [&]()->bool {
            std::streampos p = n3file.tellg();
            unsigned char b[4]; if (!n3file.read((char*)b, 4)) { n3file.clear(); n3file.seekg(p); return false; }
            n3file.clear(); n3file.seekg(p);
            return isAllowed(b[0]) && isAllowed(b[1]) && isAllowed(b[2]) && isAllowed(b[3]);
            };
        auto trySkip = [&](size_t bytes)->bool {
            std::streampos p = n3file.tellg();
            n3file.seekg(bytes, std::ios::cur);
            if (!n3file.good()) { n3file.clear(); n3file.seekg(p); return false; }
            if (peekIsFourCC()) return true;
            n3file.clear(); n3file.seekg(p); return false;
            };
        auto skipEnvelope = [&]() { for (int i = 0; i < 8; i++) (void)read_n3_value<float>(); (void)read_n3_value<int32_t>(); return true; };

        auto looks_like_fourcc = [&](uint32_t val) {
            auto b = reinterpret_cast<const unsigned char*>(&val);
            for (int i = 0; i < 4; i++) if (b[i] < 'A' || b[i] > 'Z') return false;
            return true;
            };

        auto skip_floats = [&](int n) { for (int i = 0; i < n; i++) (void)read_n3_value<float>(); };
        auto skip_ints = [&](int n) { for (int i = 0; i < n; i++) (void)read_n3_value<int32_t>(); };



        if (tag_4cc == "PGRI") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "CASH") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SHDR") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "NSKF") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SINT") { (void)read_n3_string(); (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SFLT") { (void)read_n3_string(); (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SVEC") { (void)read_n3_string(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SFRG") { int32_t g = read_n3_value<int32_t>(); int32_t n = read_n3_value<int32_t>(); for (int i = 0; i < n; i++) (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "ANIM" || tag_4cc == "VART") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "NJNT" || tag_4cc == "NJMS" || tag_4cc == "NSKL") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "JONT") { (void)read_n3_value<int32_t>(); (void)read_n3_value<int32_t>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); (void)read_n3_string(); return true; }
        else if (tag_4cc == "JOMS") { (void)read_n3_string(); { int32_t n = read_n3_value<int32_t>(); for (int i = 0; i < n; i++) (void)read_n3_value<float>(); } return true; }
        else if (tag_4cc == "SKNL") { (void)read_n3_string(); { int32_t n = read_n3_value<int32_t>(); for (int i = 0; i < n; i++) (void)read_n3_string(); } if (n3version == 2) { (void)read_n3_value<uint8_t>(); (void)read_n3_string(); } return true; }
        else if (tag_4cc == "LBOX") { for (int i = 0; i < 8; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "MNTP") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SSTA") { (void)read_n3_string(); (void)read_n3_string(); return true; }
        else if (tag_4cc == "HRCH" || tag_4cc == "HCRH") {
            (void)read_n3_value<uint8_t>();
            return true;
        }

        // Rotation (quat)
        else if (tag_4cc == "ROTN") {
            for (int i = 0; i < 4; i++) (void)read_n3_value<float>();
            return true;
        }

        // Position (point)
        else if (tag_4cc == "POSI") {
            for (int i = 0; i < 4; i++) (void)read_n3_value<float>();
            return true;
        }

        // Scale (vector)
        else if (tag_4cc == "SCAL") {
            for (int i = 0; i < 4; i++) (void)read_n3_value<float>();
            return true;
        }

        else if (tag_4cc == "RPIV" || tag_4cc == "VIPR") {
            for (int i = 0; i < 4; i++) (void)read_n3_value<float>();
            return true;
        }

        else if (tag_4cc == "SPIV" || tag_4cc == "VIPS") {
            for (int i = 0; i < 4; i++) (void)read_n3_value<float>();
            return true;
        }

        else if (tag_4cc == "SVSP" || tag_4cc == "PSVS") {
            (void)read_n3_value<uint8_t>();
            return true;
        }
        else if (tag_4cc == "SLKV" || tag_4cc == "VKLS") {
            (void)read_n3_value<uint8_t>();
            return true;
        }      

        else if (tag_4cc == "BASE" || tag_4cc == "ESAB") {
            (void)read_n3_value<int32_t>();
            return true;
        }

        else if (tag_4cc == "SLPT" || tag_4cc == "TPLS") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "ANNO" || tag_4cc == "ONNA") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SPNM" || tag_4cc == "MNPS") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SVCN" || tag_4cc == "NCVS") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SANI" || tag_4cc == "INAS") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SAGR" || tag_4cc == "RGAS") { (void)read_n3_value<int32_t>(); return true; }

        else if (tag_4cc == "ADDK" || tag_4cc == "KDDA") {
            std::string valueType = read_n3_string();   
            int32_t count = read_n3_value<int32_t>();
            if (valueType == "Float") {
                for (int i = 0; i < count; i++) { (void)read_n3_value<float>(); (void)0; }
            }
            else if (valueType == "Float4") {
                for (int i = 0; i < count; i++) { for (int j = 0; j < 4; j++) (void)read_n3_value<float>(); }
            }
            else if (valueType == "Int") {
                for (int i = 0; i < count; i++) { (void)read_n3_value<int32_t>(); }
            }
            else {
                // unknown type? do nothing
            }
            return true;
        }

        else if (tag_4cc == "ADPK" || tag_4cc == "KPDA" ||
            tag_4cc == "ADEK" || tag_4cc == "KEDA" ||
            tag_4cc == "ADSK" || tag_4cc == "KSDA") {
            int32_t count = read_n3_value<int32_t>();
            std::streampos afterCount = n3file.tellg();
            auto try_stride = [&](bool uv)->bool {
                n3file.seekg(afterCount);
                for (int i = 0; i < count; i++) {
                    if (uv) { (void)read_n3_value<int32_t>(); skip_floats(5); }   // layer + time,u,v,skip,skip
                    else { skip_floats(5); }                                    // time,x,y,z,w
                }
                std::streampos p = n3file.tellg();
                uint32_t probe = read_n3_value<uint32_t>();
                n3file.seekg(p);
                return looks_like_fourcc(probe);
                };

            if (try_stride(false)) {          // transform animator
                n3file.seekg(afterCount);
                for (int i = 0; i < count; i++) skip_floats(5);
                return true;
            }
            if (try_stride(true)) {           // UV animator
                n3file.seekg(afterCount);
                for (int i = 0; i < count; i++) { (void)read_n3_value<int32_t>(); skip_floats(5); }
                return true;
            }
            return false;
        }

        // LOD distances
        else if (tag_4cc == "SMID" || tag_4cc == "DIMS") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SMAD" || tag_4cc == "DAMS") { (void)read_n3_value<float>(); return true; }

        else if (tag_4cc == "SLCB") { for (int i = 0; i < 6; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SPOS") { for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SQUT") { for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SSCL") { for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SRTP") { for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SSCP") { for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SVSP") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SLKV") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SMID") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SMAD") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SSHD") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SMSH") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SGRI") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SANM") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "BJNT") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SJNT") { (void)read_n3_value<int32_t>(); (void)read_n3_value<int32_t>(); for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); for (int i = 0; i < 3; i++) (void)read_n3_value<float>(); (void)read_n3_string(); return true; }
        else if (tag_4cc == "SVRT") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "BGFR") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SFGI") { (void)read_n3_value<int32_t>(); (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "BGJP") { (void)read_n3_value<int32_t>(); (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SJID") { (void)read_n3_value<int32_t>(); (void)read_n3_value<int32_t>(); for (int i = 0; i < 8; i++) (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "EDFR") { return true; }

        else if (tag_4cc == "EFRQ") return skipEnvelope();
        else if (tag_4cc == "PLFT") return skipEnvelope();
        else if (tag_4cc == "PSMN") return skipEnvelope();
        else if (tag_4cc == "PSMX") return skipEnvelope();
        else if (tag_4cc == "PSVL") return skipEnvelope();
        else if (tag_4cc == "PRVL") return skipEnvelope();
        else if (tag_4cc == "PSZE") return skipEnvelope();
        else if (tag_4cc == "PMSS") return skipEnvelope();
        else if (tag_4cc == "PTMN") return skipEnvelope();
        else if (tag_4cc == "PVLF") return skipEnvelope();
        else if (tag_4cc == "PAIR") return skipEnvelope();
        else if (tag_4cc == "PRED") return skipEnvelope();
        else if (tag_4cc == "PGRN") return skipEnvelope();
        else if (tag_4cc == "PBLU") return skipEnvelope();
        else if (tag_4cc == "PALP") return skipEnvelope();

        else if (tag_4cc == "SEMD") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SLOP") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SACD") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SROF") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SBBO") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SRMN") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SRMX") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SGRV") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SPST") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "STTX") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SSTS") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SCVR") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SCVS") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SCVT") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SPCT") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SCVU") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "SSDT") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SVAF") { (void)read_n3_value<uint8_t>(); return true; }
        else if (tag_4cc == "STDL") { (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "STUS") { (void)read_n3_value<int32_t>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }
        else if (tag_4cc == "SSPI") { (void)read_n3_value<int32_t>(); for (int i = 0; i < 4; i++) (void)read_n3_value<float>(); return true; }

        else if (tag_4cc == "SCVA") return skipEnvelope();
        else if (tag_4cc == "SCVB") return skipEnvelope();
        else if (tag_4cc == "SCVD") return skipEnvelope();
        else if (tag_4cc == "SCVE") return skipEnvelope();
        else if (tag_4cc == "SCVF") return skipEnvelope();
        else if (tag_4cc == "SCVH") return skipEnvelope();
        else if (tag_4cc == "SCVJ") return skipEnvelope();
        else if (tag_4cc == "SCVL") return skipEnvelope();
        else if (tag_4cc == "STMM") return skipEnvelope();
        else if (tag_4cc == "SCVN") return skipEnvelope();
        else if (tag_4cc == "SCVQ") return skipEnvelope();
        else if (tag_4cc == "SCVC") return skipEnvelope();
        else if (tag_4cc == "SCVM") return skipEnvelope();

        else if (tag_4cc == "ADDA") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SLPT") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SPNM") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SVCN") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SANI") { (void)read_n3_string(); return true; }
        else if (tag_4cc == "SAGR") { (void)read_n3_value<int32_t>(); return true; }
        else if (tag_4cc == "SSTA") { (void)read_n3_string(); (void)read_n3_string(); return true; }
		
        return false;
    }


    bool parse_tag_state(const std::string& tag_4cc, std::shared_ptr<Node>& node) {
        if (tag_4cc == "MNMT") {
            read_n3_string();
        }
        else if (tag_4cc == "MATE") {
            read_n3_string();
        }
        else if (tag_4cc == "STXT") {
            std::string tex_type = read_n3_string();
            std::string tex_name = read_n3_string();
            node->shader_textures[tex_type] = tex_name;

            std::string full = resolveTexture(tex_name);
            GLuint texId = LoadTexture2D(full);
            if (texId) {
                std::cout << "[INFO] Bound texture " << tex_type << " -> " << texId << " (" << full << ")\n";
            }
        }
        else {
            return false;
        }
        return true;
    }

    bool parse_file(const std::string& filepath) {
        n3file.open(filepath, std::ios::binary);
        if (!n3file.is_open()) {
            report("ERROR", "Could not open file: " + filepath);
            return false;
        }

        // Header
        char header_bytes[4];
        n3file.read(header_bytes, 4);
        std::string header_4cc(header_bytes, 4);

        if (header_4cc == "NEB3") {
            byteorder = "little";
        }
        else if (header_4cc == "3BEN") {
            byteorder = "big";
        }
        else {
            report("ERROR", "Invalid file, unknown fourCC '" + header_4cc + "'");
            n3file.close();
            return false;
        }

        // Version
        uint32_t header_version;
        n3file.read(reinterpret_cast<char*>(&header_version), sizeof(header_version));
        std::cout << "n3 Version: " << header_version << std::endl;

        if (std::find(SUPPORTED_VERSIONS.begin(), SUPPORTED_VERSIONS.end(), header_version) == SUPPORTED_VERSIONS.end()) {
            report("ERROR", "Unsupported version '" + std::to_string(header_version) + "'");
            n3file.close();
            return false;
        }

        n3version = header_version;

        bool done = false;
        std::shared_ptr<Node> current_node = nullptr;
        int current_node_idx = -1;

        while (!done) {
            std::string tag_4cc = read_n3_fourcc();

            if (tag_4cc == ">MDL") {
                n3modeltype = read_n3_fourcc();
                n3modelname = read_n3_string();
                std::cout << ">MDL\n";
                std::cout << "model_type_4cc: '" << n3modeltype << "'\n";
                std::cout << "model_name: '" << n3modelname << "'\n";
            }
            else if (tag_4cc == "<MDL") {
                done = true;
                current_node = nullptr;
            }
            else if (tag_4cc == ">MND") {
                std::string node_type_4cc = read_n3_fourcc();
                std::string node_name = read_n3_string();

                auto new_node = std::make_shared<Node>(Node{ node_name, node_type_4cc, current_node });
                std::cout << ">MND " << new_node->node_type << " - " << new_node->node_name << "\n";

                if (current_node) {
                    current_node->node_children.push_back(new_node);
                }

                n3node_list.push_back(new_node);
                current_node = new_node;
                current_node_idx++;
            }
            else if (tag_4cc == "<MND") {
                current_node_idx--;
                if (current_node_idx >= 0) {
                    current_node = n3node_list[current_node_idx];
                }
                else {
                    current_node = nullptr;
                }
            }
            else if (tag_4cc == "EOF_") {
                done = true;
                current_node = nullptr;
            }
            else {
                if (current_node) {
                    if (!parse_tag_shape(tag_4cc, current_node) &&
                        !parse_tag_state(tag_4cc, current_node) &&
                        !skip4cc(tag_4cc, current_node))
                    {
                        report("ERROR", "Unknown tag '" + tag_4cc + "'");
                        n3file.close();
                        return false;
                    }
                }
                else {
                    report("ERROR", "Tag without active node: " + tag_4cc);
                    n3file.close();
                    return false;
                }
            }
        }

        n3file.close();
        return true;
    }
};

template<typename T>
T clamp(T v, T lo, T hi) {
    return (v < lo) ? lo : (v > hi) ? hi : v;
}
 
GLuint gridVAO = 0, gridVBO = 0;
GLsizei gridVertexCount = 0;
constexpr double M_PI = 3.14159265358979323846;

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
        pitch -= (float)dy * 0.01f;  // Inverted for natural mouse movement
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
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    distance -= yoffset * 0.5f;
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


struct Vec3 {
    float x, y, z;
    Vec3() : x(0), y(0), z(0) {}
    Vec3(float x_, float y_, float z_) : x(x_), y(y_), z(z_) {}
};
struct ObjVertex {
    float px, py, pz;
    float nx, ny, nz;
    float tx, ty, tz;
    float bx, by, bz;
    float u0, v0;
    float u1, v1;
    float cr, cg, cb, ca;
};

struct Mesh {
    std::vector<ObjVertex> verts;
    std::vector<uint32_t> idx;
    GLuint vao = 0, vbo = 0, ebo = 0;
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


static bool LoadNVX2_AndNormalizeUnit(const std::string& path, Mesh& out) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        std::cerr << "[ERROR] Could not open NVX2 file: " << path << std::endl;
        char buf[512];
        if (_getcwd(buf, sizeof(buf))) std::cerr << "[INFO] Current working directory: " << buf << std::endl;
        return false;
    }

    char magic[4];
    f.read(magic, 4);
    if (strncmp(magic, "NVX2", 4) != 0 && strncmp(magic, "2XVN", 4) != 0) {
        std::cerr << "[ERROR] Invalid NVX2 magic header: " << magic << std::endl;
        return false;
    }

    uint32_t numGroups, numVertices, vertexWidth, numTriangles, numEdges, compMask;
    f.read(reinterpret_cast<char*>(&numGroups), 4);
    f.read(reinterpret_cast<char*>(&numVertices), 4);
    f.read(reinterpret_cast<char*>(&vertexWidth), 4);
    f.read(reinterpret_cast<char*>(&numTriangles), 4);
    f.read(reinterpret_cast<char*>(&numEdges), 4);
    f.read(reinterpret_cast<char*>(&compMask), 4);

    std::cout << "[INFO] NVX2 header: vertices=" << numVertices
        << " tris=" << numTriangles
        << " stride=" << vertexWidth * 4
        << " compMask=0x" << std::hex << compMask << std::dec << "\n";

    if (numGroups > 0) f.seekg(numGroups * 24, std::ios::cur);

    struct ParsedComponent { uint32_t type; int offset; };
    std::vector<ParsedComponent> components;
    int currentOffset = 0;

    for (int bit = 0; bit < 32; bit++) {
        uint32_t flag = 1u << bit;
        if (compMask & flag) {
            components.push_back({ flag, currentOffset });
            currentOffset += GetComponentSize(flag);
        }
    }

    int posOffset = -1, uv0Offset = -1;
    bool uv0Short = false;

    int nrmOffset = -1, tanOffset = -1, binOffset = -1, colOffset = -1, uv1Offset = -1;
    bool nrmUB4 = false, tanUB4 = false, binUB4 = false, colUB4 = false, uv1Short = false;

    for (const auto& c : components) {
        if (c.type & Coord) { posOffset = c.offset; }
        else if (c.type & Uv0) { uv0Offset = c.offset; uv0Short = false; }
        else if (c.type & Uv0S2) { uv0Offset = c.offset; uv0Short = true; }
        else if (c.type & Normal) { nrmOffset = c.offset; nrmUB4 = false; }
        else if (c.type & NormalUB4N) { nrmOffset = c.offset; nrmUB4 = true; }
        else if (c.type & Tangent) { tanOffset = c.offset; tanUB4 = false; }
        else if (c.type & TangentUB4N) { tanOffset = c.offset; tanUB4 = true; }
        else if (c.type & Binormal) { binOffset = c.offset; binUB4 = false; }
        else if (c.type & BinormalUB4N) { binOffset = c.offset; binUB4 = true; }
        else if (c.type & Color) { colOffset = c.offset; colUB4 = false; }
        else if (c.type & ColorUB4N) { colOffset = c.offset; colUB4 = true; }
        else if (c.type & Uv1) { uv1Offset = c.offset; uv1Short = false; }
        else if (c.type & Uv1S2) { uv1Offset = c.offset; uv1Short = true; }
    }

    if (posOffset < 0) {
        std::cerr << "[ERROR] No Coord in compMask\n";
        return false;
    }

    if (currentOffset != static_cast<int>(vertexWidth * 4)) {
        std::cerr << "[WARNING] Mismatched strides - calculated=(" << currentOffset
            << ") vs file header=(" << vertexWidth * 4 << ")\n";
    }

    std::vector<ObjVertex> verts(numVertices);
    std::vector<uint8_t> vbuf(vertexWidth * 4);

    auto ub4_to_n11 = [](uint8_t b) { return (float)b / 127.5f - 1.0f; };
    auto ub4_to_01 = [](uint8_t b) { return (float)b / 255.0f; };

    for (uint32_t i = 0; i < numVertices; i++) {
        f.read(reinterpret_cast<char*>(vbuf.data()), vertexWidth * 4);

        const float* pos = reinterpret_cast<const float*>(vbuf.data() + posOffset);
        verts[i].px = pos[0]; verts[i].py = pos[1]; verts[i].pz = pos[2];

        if (uv0Offset >= 0) {
            if (!uv0Short) {
                const float* uv = reinterpret_cast<const float*>(vbuf.data() + uv0Offset);
                verts[i].u0 = uv[0]; verts[i].v0 = 1.0f - uv[1];
            }
            else {
                const int16_t* uv = reinterpret_cast<const int16_t*>(vbuf.data() + uv0Offset);
                verts[i].u0 = (float)uv[0] / 8191.0f;
                verts[i].v0 = 1.0f - (float)uv[1] / 8191.0f;
            }
        }
        else { verts[i].u0 = verts[i].v0 = 0.0f; }

        if (uv1Offset >= 0) {
            if (!uv1Short) {
                const float* uv = reinterpret_cast<const float*>(vbuf.data() + uv1Offset);
                verts[i].u1 = uv[0]; verts[i].v1 = 1.0f - uv[1];
            }
            else {
                const int16_t* uv = reinterpret_cast<const int16_t*>(vbuf.data() + uv1Offset);
                verts[i].u1 = (float)uv[0] / 8191.0f;
                verts[i].v1 = 1.0f - (float)uv[1] / 8191.0f;
            }
        }
        else { verts[i].u1 = verts[i].v1 = 0.0f; }

        if (nrmOffset >= 0) {
            if (!nrmUB4) {
                const float* n = reinterpret_cast<const float*>(vbuf.data() + nrmOffset);
                verts[i].nx = n[0]; verts[i].ny = n[1]; verts[i].nz = n[2];
            }
            else {
                const uint8_t* n = vbuf.data() + nrmOffset;
                verts[i].nx = ub4_to_n11(n[0]); verts[i].ny = ub4_to_n11(n[1]); verts[i].nz = ub4_to_n11(n[2]);
            }
        }
        else { verts[i].nx = 0; verts[i].ny = 0; verts[i].nz = 1; }

        if (tanOffset >= 0) {
            if (!tanUB4) {
                const float* t = reinterpret_cast<const float*>(vbuf.data() + tanOffset);
                verts[i].tx = t[0]; verts[i].ty = t[1]; verts[i].tz = t[2];
            }
            else {
                const uint8_t* t = vbuf.data() + tanOffset;
                verts[i].tx = ub4_to_n11(t[0]); verts[i].ty = ub4_to_n11(t[1]); verts[i].tz = ub4_to_n11(t[2]);
            }
        }
        else { verts[i].tx = 1; verts[i].ty = 0; verts[i].tz = 0; }

        if (binOffset >= 0) {
            if (!binUB4) {
                const float* b = reinterpret_cast<const float*>(vbuf.data() + binOffset);
                verts[i].bx = b[0]; verts[i].by = b[1]; verts[i].bz = b[2];
            }
            else {
                const uint8_t* b = vbuf.data() + binOffset;
                verts[i].bx = ub4_to_n11(b[0]); verts[i].by = ub4_to_n11(b[1]); verts[i].bz = ub4_to_n11(b[2]);
            }
        }
        else { verts[i].bx = 0; verts[i].by = 1; verts[i].bz = 0; }

        if (colOffset >= 0) {
            if (!colUB4) {
                const float* c = reinterpret_cast<const float*>(vbuf.data() + colOffset);
                verts[i].cr = c[0]; verts[i].cg = c[1]; verts[i].cb = c[2]; verts[i].ca = c[3];
            }
            else {
                const uint8_t* c = vbuf.data() + colOffset;
                verts[i].cr = ub4_to_01(c[0]); verts[i].cg = ub4_to_01(c[1]); verts[i].cb = ub4_to_01(c[2]); verts[i].ca = ub4_to_01(c[3]);
            }
        }
        else { verts[i].cr = verts[i].cg = verts[i].cb = verts[i].ca = 1.0f; }

        if (i < 5) {
            std::cout << "[DEBUG] v" << i << " pos(" << verts[i].px << "," << verts[i].py << "," << verts[i].pz
                << ") uv0(" << verts[i].u0 << "," << verts[i].v0 << ") uv1(" << verts[i].u1 << "," << verts[i].v1
                << ") n(" << verts[i].nx << "," << verts[i].ny << "," << verts[i].nz << ")\n";
        }
    }

    out.idx.resize(numTriangles * 3);
    for (uint32_t i = 0; i < numTriangles; i++) {
        uint16_t a, b, c;
        f.read(reinterpret_cast<char*>(&a), 2);
        f.read(reinterpret_cast<char*>(&b), 2);
        f.read(reinterpret_cast<char*>(&c), 2);
        out.idx[i * 3 + 0] = a;
        out.idx[i * 3 + 1] = b;
        out.idx[i * 3 + 2] = c;
    }

    out.verts = verts;

    glGenVertexArrays(1, &out.vao);
    glGenBuffers(1, &out.vbo);
    glGenBuffers(1, &out.ebo);

    glBindVertexArray(out.vao);

    glBindBuffer(GL_ARRAY_BUFFER, out.vbo);
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(ObjVertex), verts.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, out.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, out.idx.size() * sizeof(uint32_t), out.idx.data(), GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, px));

    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, u0));

    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, nx));

    glEnableVertexAttribArray(3);
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, tx));

    glEnableVertexAttribArray(4);
    glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, bx));

    glEnableVertexAttribArray(5);
    glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, cr));

    glEnableVertexAttribArray(6);
    glVertexAttribPointer(6, 2, GL_FLOAT, GL_FALSE, sizeof(ObjVertex), (void*)offsetof(ObjVertex, u1));

    glBindVertexArray(0);
    return true;
}

// Vertex shader
const char* vertexShaderSource = R"(
#version 330 core
layout (location=0) in vec3 aPos;
layout (location=1) in vec2 aUV0;
layout (location=2) in vec3 aNormal;
layout (location=3) in vec3 aTangent;
// layout (location=4) in vec3 aBinormal;

out VS_OUT {
    vec2 uv;
    mat3 TBN;
    vec3 wPos;
    vec3 wN;
} vs;

uniform mat4 model, view, projection;

void main() {
    mat3 m3 = mat3(model);
    vec3 N = normalize(m3 * aNormal);
    vec3 T = normalize(m3 * aTangent);
    vec3 B = normalize(cross(N, T));

    vs.TBN = mat3(T, B, N);
    vs.uv = aUV0;
    vec4 wp = model * vec4(aPos, 1.0);
    vs.wPos = wp.xyz;
    vs.wN = N;
    gl_Position = projection * view * wp;
}

)";

// Fragment shader
const char* fragmentShaderSource = R"(
  #version 330 core
in VS_OUT {
    vec2 uv;
    mat3 TBN;
    vec3 wPos;
    vec3 wN;
} fs;

out vec4 FragColor;

uniform sampler2D DiffuseMap;
uniform sampler2D SpecMap;
uniform sampler2D BumpMap;
uniform sampler2D EmissiveMap;

uniform vec3  CameraPos;
uniform float SpecularPower;
uniform float SpecularIntensity;
uniform float EmissiveIntensity;
uniform bool  UseBump;

vec3 decodeNormalAuto(vec4 nm, mat3 TBN, vec3 wN) {
    vec2 ag = nm.ag * 2.0 - 1.0;                  // DXT5nm (A/G)
    vec2 rg = nm.rg * 2.0 - 1.0;                  // BC5/ATI2 (R/G)
    float zag = sqrt(max(1.0 - dot(ag, ag), 0.0));
    float zrg = sqrt(max(1.0 - dot(rg, rg), 0.0));
    vec3 N_ag = normalize(TBN * vec3(ag, zag));
    vec3 N_rg = normalize(TBN * vec3(rg, zrg));
    return (dot(N_ag, wN) >= dot(N_rg, wN)) ? N_ag : N_rg;
}

void main() {
    vec3 lightDir = normalize(vec3(-0.5, -1.0, -0.3));

    vec3 albedo = texture(DiffuseMap, fs.uv).rgb;

    vec3 N = fs.wN;
    if (UseBump) {
        vec4 nm = texture(BumpMap, fs.uv);
        N = decodeNormalAuto(nm, fs.TBN, fs.wN);
    }

    vec3 L = normalize(-lightDir);
    vec3 V = normalize(CameraPos - fs.wPos);
    vec3 H = normalize(L + V);

    float NdotL = max(dot(N, L), 0.0);
    float specMask = texture(SpecMap, fs.uv).r;
    float spec = pow(max(dot(N, H), 0.0), SpecularPower) * SpecularIntensity * specMask;

    vec3 emissive = texture(EmissiveMap, fs.uv).rgb * EmissiveIntensity;
    vec3 ambient = 1.0 * albedo;
    vec3 color = ambient + albedo * NdotL + vec3(spec) + emissive;
    FragColor = vec4(color, 1.0);
}

)";

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

void setIdentity(float* matrix) {
    for (int i = 0; i < 16; i++) matrix[i] = 0.0f;
    matrix[0] = matrix[5] = matrix[10] = matrix[15] = 1.0f;
}

void multiplyMatrices(float* result, const float* a, const float* b) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i * 4 + j] = 0.0f;
            for (int k = 0; k < 4; k++) {
                result[i * 4 + j] += a[i * 4 + k] * b[k * 4 + j];
            }
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
    float f = 1.0f / tan(fov * M_PI / 180.0f / 2.0f);
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

int main(int argc, char* argv[])
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    double sizeX = 1200, sizeY = 800;
    GLFWwindow* window = glfwCreateWindow(sizeX, sizeY, "Simo's Neb3 Viewer", NULL, NULL);
    if (!window) { glfwTerminate(); return -1; }

    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetKeyCallback(window, [](GLFWwindow* w, int key, int sc, int action, int mods) {
        if (key >= 0 && key < 1024) { if (action == GLFW_PRESS) keys[key] = true; else if (action == GLFW_RELEASE) keys[key] = false; }
        });
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

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

    glGenVertexArrays(1, &gridVAO);
    glGenBuffers(1, &gridVBO);
    glBindVertexArray(gridVAO);
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(gridVertices), gridVertices, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

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

    InitGrid(10.0f, 1.0f);

    std::vector<std::string> n3Files = {
        "C:/drasa_online/work/models/arch_house_01.n3",
        "C:/drasa_online/work/models/arch_well_01.n3",
        "C:/drasa_online/work/models/apple_fresh_green.n3"
    };

    struct Loaded {
        Mesh mesh;
        GLuint texDiffuse = 0, texSpec = 0, texBump = 0, texEmissive = 0;
        float pos[3];
    };

    std::vector<Loaded> sceneObjects;
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<float> dist(-10.0f, 10.0f);

    for (auto& n3FullPath : n3Files) {
        Parser parser;
        if (!parser.parse_file(n3FullPath)) continue;

        auto resolveDDS = [](std::string t)->std::string {
            if (t.rfind("tex:", 0) == 0) t = t.substr(4);
            return "C:/drasa_online/work/textures/" + t + ".dds";
            };
        auto loadTexByKey = [&](const char* key)->GLuint {
            for (auto& n : parser.getNodes()) {
                auto it = n->shader_textures.find(key);
                if (it != n->shader_textures.end())
                    return LoadTexture2D(resolveDDS(it->second));
            }
            return 0u;
            };

        Loaded obj;
        obj.texDiffuse = loadTexByKey("DiffMap0");
        obj.texSpec = loadTexByKey("SpecMap0");
        obj.texBump = loadTexByKey("BumpMap0");
        obj.texEmissive = loadTexByKey("EmsvMap0");
        obj.pos[0] = dist(rng);
        obj.pos[1] = 0.0f;
        obj.pos[2] = dist(rng);
        sceneObjects.push_back(std::move(obj));

        std::string nvx2Path;
        for (auto& node : parser.getNodes()) {
            if (!node->mesh_resource_id.empty()) {
                nvx2Path = node->mesh_resource_id;
                if (nvx2Path.rfind("msh:", 0) == 0) nvx2Path = nvx2Path.substr(4);
                nvx2Path = "C:/drasa_online/work/meshes/" + nvx2Path;
                break;
            }
        }
        if (nvx2Path.empty()) continue;
        if (!LoadNVX2_AndNormalizeUnit(nvx2Path, obj.mesh)) continue;

        sceneObjects.push_back(std::move(obj));
    }

    glUseProgram(shaderProgram);
    glUniform1i(glGetUniformLocation(shaderProgram, "DiffuseMap"), 0);
    glUniform1i(glGetUniformLocation(shaderProgram, "SpecMap"), 1);
    glUniform1i(glGetUniformLocation(shaderProgram, "BumpMap"), 2);
    glUniform1i(glGetUniformLocation(shaderProgram, "EmissiveMap"), 3);
    glUniform1f(glGetUniformLocation(shaderProgram, "SpecularPower"), 26.0f);
    glUniform1f(glGetUniformLocation(shaderProgram, "SpecularIntensity"), 1.0f);
    glUniform1f(glGetUniformLocation(shaderProgram, "EmissiveIntensity"), 1.0f);
    glUniform1i(glGetUniformLocation(shaderProgram, "UseBump"), GL_TRUE);
    glUseProgram(0);

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    while (!glfwWindowShouldClose(window))
    {
        processInput(window);

        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glUseProgram(shaderProgram);

        float projectionMatrix[16];
        createPerspectiveMatrix(projectionMatrix, 45.0f, (float)(sizeX / sizeY), 0.1f, 100.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1, GL_FALSE, projectionMatrix);

        float viewMatrix[16];
        float camX = sinf(yaw) * cosf(pitch) * distance + targetX;
        float camY = sinf(pitch) * distance + targetY;
        float camZ = cosf(yaw) * cosf(pitch) * distance + targetZ;
        createLookAtMatrix(viewMatrix, camX, camY, camZ, targetX, targetY, targetZ, 0.0f, 1.0f, 0.0f);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, GL_FALSE, viewMatrix);
        glUniform3f(glGetUniformLocation(shaderProgram, "CameraPos"), camX, camY, camZ);

        float modelMatrix[16];
        setIdentity(modelMatrix);
        glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, modelMatrix);

        for (auto& obj : sceneObjects) {
            float modelMatrix[16];
            setIdentity(modelMatrix);
            translateMatrix(modelMatrix, obj.pos[0], obj.pos[1], obj.pos[2]);
            glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "model"), 1, GL_FALSE, modelMatrix);

            glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, obj.texDiffuse);
            glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, obj.texSpec);
            glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, obj.texBump);
            glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, obj.texEmissive);

            glBindVertexArray(obj.mesh.vao);
            glDrawElements(GL_TRIANGLES, (GLsizei)obj.mesh.idx.size(), GL_UNSIGNED_INT, 0);
            glBindVertexArray(0);
        }

        DrawGrid(gridShaderProgram, viewMatrix, projectionMatrix);

        glUseProgram(0);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    CleanupGrid();
    for (auto& obj : sceneObjects) {
        glDeleteVertexArrays(1, &obj.mesh.vao);
        glDeleteBuffers(1, &obj.mesh.vbo);
        glDeleteBuffers(1, &obj.mesh.ebo);
        if (obj.texDiffuse)  glDeleteTextures(1, &obj.texDiffuse);
        if (obj.texSpec)     glDeleteTextures(1, &obj.texSpec);
        if (obj.texBump)     glDeleteTextures(1, &obj.texBump);
        if (obj.texEmissive) glDeleteTextures(1, &obj.texEmissive);
    }
    sceneObjects.clear();
    glDeleteProgram(shaderProgram);
    glDeleteProgram(gridShaderProgram);
    glfwTerminate();
    return 0;
}

