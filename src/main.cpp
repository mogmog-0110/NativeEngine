#include <cstdio>
#include <cstring>

namespace ne { int runSelftest(); }

int main(int argc, char** argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) {
        return ne::runSelftest();
    }
    std::printf("NativeEngine -- native-PBC rigid-body engine for macrocycle assembly\n");
    std::printf("usage:\n");
    std::printf("  NativeEngine --selftest    run the physics unit tests\n");
    return 0;
}
