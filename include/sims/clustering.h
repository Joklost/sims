#ifndef MANETSIMS_CLUSTERING_H
#define MANETSIMS_CLUSTERING_H

#include <vector>
#include <limits>
#include <cmath>
#include <queue>
#include <unordered_map>

#include "node.h"

namespace sims {
    constexpr double UNDEFINED = -0.1;

    using neighbour_t = ::std::pair<Node &, double>;
    using neighbourhood_t = ::std::vector<neighbour_t>;

    class Optics {
    private:
        struct Neighbour {
            uint32_t node{};
            double distance{};
        };

    public:
        class Cluster {
        public:
            Cluster(uint32_t id, ::std::vector<Node> &nodes) : id(id), nodes(nodes) {}

            geo::Location centroid() const;

            unsigned long size() const;

            double radius() const;

            uint32_t get_id() const;

            double cost() const;

            const ::std::vector<Node> &get_nodes() const;

            bool contains(const Node &node) const;

            bool operator==(const Cluster &rhs) const;

            bool operator!=(const Cluster &rhs) const;

        private:
            uint32_t id{};
            ::std::vector<Node> nodes{};

            mutable bool cached = false;
            mutable geo::Location _centroid{};
        };

        class CLink {
        public:
            CLink(uint64_t id, Cluster &c1, Cluster &c2);

            const ::std::pair<Cluster, Cluster> &get_clusters() const;

            double get_distance() const;

            uint64_t get_id() const;

            bool operator==(const CLink &rhs) const;

            bool operator!=(const CLink &rhs) const;

            bool operator<(const CLink &rhs) const;

            bool operator>(const CLink &rhs) const;

            bool operator<=(const CLink &rhs) const;

            bool operator>=(const CLink &rhs) const;

            bool contains(const Node &node) const;

            bool contains(const Node &n1, const Node &n2) const;

            double get_rssi() const;

            void set_rssi(double rssi);

            double get_pep() const;

            void set_pep(double pep);

        private:
            double distance;
            double rssi;
            double pep;
            uint64_t id;
            ::std::pair<Cluster, Cluster> clusters;

        };

        Optics();

        ::std::vector<Node> compute_ordering(::std::vector<Node> &nodes, double eps /* kilometers */, int minpts);

        ::std::vector<Cluster> cluster(::std::vector<Node> &ordering);

        ::std::vector<Cluster> cluster(::std::vector<Node> &ordering, double threshold);



private:
        void clusterize_remaining(const std::vector<Node> &nodes, std::vector<Optics::Cluster> &clusters);

        void processed(Node &p);

        double core_distance(Node &p);

        ::std::vector<Neighbour> &compute_neighbours(Node &p);

        void update_seeds(Node &p, ::std::vector<uint32_t> &seeds);

        double eps;
        int minpts;

        ::std::unordered_map<uint32_t, Node> graph{};
        ::std::vector<int> unprocessed{};
        ::std::vector<Node> ordered{};
        ::std::unordered_map<Node, ::std::vector<Neighbour>> neighbourhoods{};
    };

}


#endif /* MANETSIMS_CLUSTERING_H */
