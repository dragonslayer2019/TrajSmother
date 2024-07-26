// created by yuqing.wu on 03/06/2024

#include <string>

#ifndef TYPE_H
#define TYPE_H

namespace global_path_search {

enum class CubeType {CUBE_ONE = 0, CUBE_TWO = 1};
inline std::string ToString(CubeType cube_type) {
    switch (cube_type) {
        case CubeType::CUBE_ONE:
            return "Cube-One";
        case CubeType::CUBE_TWO:
            return "Cube-Two";
    }
    return "Unknown Cube";
}

enum class VerticeType {CENTER = 0, SURFACE = 1};
inline std::string ToString(VerticeType vertice_type) {
    switch (vertice_type) {
        case VerticeType::CENTER:
            return "center";
        case VerticeType::SURFACE:
            return "surface";
    }
    return "Unknown vertice";
}

enum class SearchResult { SUCCESS = 0, FAILURE = 1, RUNNING = 2 };
inline std::string ToString(SearchResult res) {
    switch (res) {
        case SearchResult::SUCCESS:
            return "success";
        case SearchResult::FAILURE:
            return "failure";
        case SearchResult::RUNNING:
            return "running";
    }
    return "unknown result";
}

enum class NodeState { IN_OPEN_SET = 0, IN_CLOSET_SET = 1, UNVISITED = 2 };
inline std::string ToString(NodeState state) {
    switch (state) {
        case NodeState::IN_OPEN_SET:
            return "in open set";
        case NodeState::IN_CLOSET_SET:
            return "in close set";
        case NodeState::UNVISITED:
            return "unvisited";
    }
    return "unknown state";
}

enum class OctNodeType { SRC = 0, DES = 1, NORMAL = 2 };
inline std::string ToString(OctNodeType node_type) {
    switch (node_type) {
        case OctNodeType::SRC:
            return "source";
        case OctNodeType::DES:
            return "destination";
        case OctNodeType::NORMAL:
            return "normal";
    }
    return "unknown octnode type";
}

} // namespace global_path_search

#endif // TYPE_H