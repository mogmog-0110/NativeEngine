#include <cstdio>
#include <cstring>

namespace ne {
int runSelftest();
}

int main(int argc, char** argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) {
        return ne::runSelftest();
    }
    std::printf("NativeEngine -- deterministic rigid-body engine with native PBC\n");
    std::printf("usage:\n");
    std::printf("  NativeEngine --selftest    run the physics unit tests\n");
    std::printf("\nInteractive demos live in the viewer (build_viewer.bat):\n");
    std::printf("  NativeViewer <1..9>        live, interactive PhysX-style demo scenes\n");
    return 0;
}
