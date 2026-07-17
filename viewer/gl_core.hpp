#pragma once
#ifndef NATIVEENGINE_VIEWER_GL_CORE_HPP
#define NATIVEENGINE_VIEWER_GL_CORE_HPP

// Minimal OpenGL 3.3-core loader on top of freeglut. freeglut ships only the GL
// 1.1 headers, so the modern entry points (VAOs, VBOs, shaders, instancing) are
// loaded at runtime via glutGetProcAddress. We declare exactly the functions the
// viewer uses, plus the few enum values, to avoid needing glext.h.
#include <GL/freeglut.h>
#include <cstddef>

#ifndef APIENTRY
#define APIENTRY __stdcall
#endif

// -- enums (standard GL values) --
#define GL_ARRAY_BUFFER 0x8892
#define GL_ELEMENT_ARRAY_BUFFER 0x8893
#define GL_STATIC_DRAW 0x88E4
#define GL_DYNAMIC_DRAW 0x88E8
#define GL_FRAGMENT_SHADER 0x8B30
#define GL_VERTEX_SHADER 0x8B31
#define GL_COMPILE_STATUS 0x8B81
#define GL_LINK_STATUS 0x8B82
#define GL_MULTISAMPLE 0x809D
#define GL_FRAMEBUFFER_SRGB 0x8DB9

typedef char GLchar;
typedef ptrdiff_t GLsizeiptr;
typedef ptrdiff_t GLintptr;

// -- function pointer typedefs + externs --
typedef GLuint(APIENTRY* PFN_glCreateShader)(GLenum);
typedef void(APIENTRY* PFN_glShaderSource)(GLuint, GLsizei, const GLchar* const*, const GLint*);
typedef void(APIENTRY* PFN_glCompileShader)(GLuint);
typedef void(APIENTRY* PFN_glGetShaderiv)(GLuint, GLenum, GLint*);
typedef void(APIENTRY* PFN_glGetShaderInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
typedef void(APIENTRY* PFN_glDeleteShader)(GLuint);
typedef GLuint(APIENTRY* PFN_glCreateProgram)(void);
typedef void(APIENTRY* PFN_glAttachShader)(GLuint, GLuint);
typedef void(APIENTRY* PFN_glLinkProgram)(GLuint);
typedef void(APIENTRY* PFN_glGetProgramiv)(GLuint, GLenum, GLint*);
typedef void(APIENTRY* PFN_glGetProgramInfoLog)(GLuint, GLsizei, GLsizei*, GLchar*);
typedef void(APIENTRY* PFN_glUseProgram)(GLuint);
typedef void(APIENTRY* PFN_glGenVertexArrays)(GLsizei, GLuint*);
typedef void(APIENTRY* PFN_glBindVertexArray)(GLuint);
typedef void(APIENTRY* PFN_glGenBuffers)(GLsizei, GLuint*);
typedef void(APIENTRY* PFN_glBindBuffer)(GLenum, GLuint);
typedef void(APIENTRY* PFN_glBufferData)(GLenum, GLsizeiptr, const void*, GLenum);
typedef void(APIENTRY* PFN_glBufferSubData)(GLenum, GLintptr, GLsizeiptr, const void*);
typedef void(APIENTRY* PFN_glVertexAttribPointer)(GLuint, GLint, GLenum, GLboolean, GLsizei, const void*);
typedef void(APIENTRY* PFN_glEnableVertexAttribArray)(GLuint);
typedef void(APIENTRY* PFN_glVertexAttribDivisor)(GLuint, GLuint);
typedef GLint(APIENTRY* PFN_glGetUniformLocation)(GLuint, const GLchar*);
typedef void(APIENTRY* PFN_glUniformMatrix4fv)(GLint, GLsizei, GLboolean, const GLfloat*);
typedef void(APIENTRY* PFN_glUniform3fv)(GLint, GLsizei, const GLfloat*);
typedef void(APIENTRY* PFN_glUniform1f)(GLint, GLfloat);
typedef void(APIENTRY* PFN_glUniform1i)(GLint, GLint);
typedef void(APIENTRY* PFN_glDrawElementsInstanced)(GLenum, GLsizei, GLenum, const void*, GLsizei);

#define NV_GL_FUNCS(X) \
    X(PFN_glCreateShader, glCreateShader) \
    X(PFN_glShaderSource, glShaderSource) \
    X(PFN_glCompileShader, glCompileShader) \
    X(PFN_glGetShaderiv, glGetShaderiv) \
    X(PFN_glGetShaderInfoLog, glGetShaderInfoLog) \
    X(PFN_glDeleteShader, glDeleteShader) \
    X(PFN_glCreateProgram, glCreateProgram) \
    X(PFN_glAttachShader, glAttachShader) \
    X(PFN_glLinkProgram, glLinkProgram) \
    X(PFN_glGetProgramiv, glGetProgramiv) \
    X(PFN_glGetProgramInfoLog, glGetProgramInfoLog) \
    X(PFN_glUseProgram, glUseProgram) \
    X(PFN_glGenVertexArrays, glGenVertexArrays) \
    X(PFN_glBindVertexArray, glBindVertexArray) \
    X(PFN_glGenBuffers, glGenBuffers) \
    X(PFN_glBindBuffer, glBindBuffer) \
    X(PFN_glBufferData, glBufferData) \
    X(PFN_glBufferSubData, glBufferSubData) \
    X(PFN_glVertexAttribPointer, glVertexAttribPointer) \
    X(PFN_glEnableVertexAttribArray, glEnableVertexAttribArray) \
    X(PFN_glVertexAttribDivisor, glVertexAttribDivisor) \
    X(PFN_glGetUniformLocation, glGetUniformLocation) \
    X(PFN_glUniformMatrix4fv, glUniformMatrix4fv) \
    X(PFN_glUniform3fv, glUniform3fv) \
    X(PFN_glUniform1f, glUniform1f) \
    X(PFN_glUniform1i, glUniform1i) \
    X(PFN_glDrawElementsInstanced, glDrawElementsInstanced)

#define X(type, name) inline type name = nullptr;
NV_GL_FUNCS(X)
#undef X

inline bool loadGL() {
    bool ok = true;
#define X(type, name) \
    name = reinterpret_cast<type>(glutGetProcAddress(#name)); \
    if (!name) ok = false;
    NV_GL_FUNCS(X)
#undef X
    return ok;
}

#endif  // NATIVEENGINE_VIEWER_GL_CORE_HPP
