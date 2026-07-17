#include <cstdio>
#include <cstring>
#include <string>

namespace ne {
int runSelftest();
int runDemo(const std::string& name, const std::string& out);
int verifyPxrf(const std::string& path);
}

int main(int argc, char** argv) {
    if (argc >= 2 && std::strcmp(argv[1], "--selftest") == 0) {
        return ne::runSelftest();
    }
    if (argc >= 3 && std::strcmp(argv[1], "demo") == 0) {
        std::string name = argv[2];
        std::string out = (argc >= 4) ? argv[3] : ("demo_" + name + ".pxrf");
        return ne::runDemo(name, out);
    }
    if (argc >= 3 && std::strcmp(argv[1], "verifypxrf") == 0) {
        return ne::verifyPxrf(argv[2]);
    }
    std::printf("NativeEngine -- deterministic rigid-body engine with native PBC\n");
    std::printf("usage:\n");
    std::printf("  NativeEngine --selftest             run the physics unit tests\n");
    std::printf("  NativeEngine demo <name> [out.pxrf] run a demo, record a .pxrf\n");
    std::printf("    demos: pile cylinders pendulum gas pbc_gas pbc_pair all\n");
    return 0;
}
