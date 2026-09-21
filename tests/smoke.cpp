#include "rar/RemeshConfig.h"

#include <cassert>

int main() {
    rar::RemeshConfig cfg;
    assert(cfg.epsilon > 0.0);
    assert(cfg.min_edge_length > 0.0);
    assert(cfg.max_edge_length >= cfg.min_edge_length);
    assert(cfg.iterations > 0);
    return 0;
}
