#ifndef CZINSERTIONTOOL_H
#define CZINSERTIONTOOL_H

#include "node.h"
#include "mesh.h"
#include <tuple>
#include <algorithm>
#include <spdlog/spdlog.h>

namespace icy {class CZInsertionTool;}

class icy::CZInsertionTool
{
public:
    struct NodeExtension
    {
        std::vector<Element*> adj_elems;
        std::vector<int> adj_grains;
    };

    struct Facet
    {
        icy::Element *elems[2] {};
        int facet_idx[2] {};
        int orientation=-77; // index of the node in facet 1 that matches node 0 in facet 0
        std::tuple<int,int,int> key;
        static std::tuple<int,int,int> make_key(icy::Node *nd0, icy::Node *nd1, icy::Node *nd2)
        {
            int nds[3] = {nd0->globId,nd1->globId,nd2->globId};
            std::sort(std::begin(nds),std::end(nds));
            std::tuple<int,int,int> result(nds[0],nds[1],nds[2]);
            return result;
        }
    };



    void InsertCZs(icy::Mesh &mesh);

private:
    void ReplaceNodeInElement(icy::Element *elem, icy::Node *whichNode, icy::Node *replacement);
};

#endif // CZINSERTIONTOOL_H
