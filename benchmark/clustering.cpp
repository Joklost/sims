#include <iostream>
#include <iomanip>
#include <future>
#include <cmath>

#include <sims/node.h>
#include <sims/datagen.h>
#include <sims/math.h>
#include <sims/cholesky.h>
#include <sims/radiomodel.h>

#include <mpilib/helpers.h>
#include <mpilib/httpclient.h>


#define MIN_PTS 2
#define MAX_PTS 10
#define PTS_STEPS 1
#define MIN_EPS 0.05
#define MAX_EPS 1.01
#define MAX_EPS_STEPS 46
#define EPS_STEP 10_m
#define TX_DBM 26.0
#define AREA 1.0
#define NODES 1000

template<class...Durations, class DurationIn>
std::tuple<Durations...> break_down_durations(DurationIn d) {
    std::tuple<Durations...> retval;
    using discard=int[];
    (void) discard{0, (void(((std::get<Durations>(retval) = std::chrono::duration_cast<Durations>(d)),
            (d -= std::chrono::duration_cast<DurationIn>(std::get<Durations>(retval))))), 0)...};
    return retval;
}


std::string format_duration(std::chrono::microseconds us) {
    auto dur = break_down_durations<std::chrono::seconds, std::chrono::milliseconds, std::chrono::microseconds>(us);
    std::ostringstream oss;
    oss << std::setfill('0')
        << std::get<0>(dur).count()
        << "::"
        << std::setw(3)
        << std::get<1>(dur).count()
        << "::"
        << std::setw(3)
        << std::get<2>(dur).count();
    return oss.str();
}

/* This block is for clustering experiments on distance threshold */
int main(int argc, char *argv[]) {
    geo::Location upper{57.01266813458001, 10.994625734716218};
    auto lower = geo::square(upper, 10_km);
    auto nodes = sims::data::generate_nodes(NODES, upper, lower);

    std::vector<double> thresholds{100_m, 125_m, 150_m, 175_m, 200_m, 225_m, 250_m, 275_m, 300_m, 325_m, 350_m, 375_m,
                                   400_m, 425_m, 450_m, 475_m, 500_m, 525_m, 550_m, 575_m, 600_m, 625_m, 650_m, 675_m,
                                   700_m, 725_m, 750_m, 775_m, 800_m, 825_m, 850_m, 875_m, 900_m, 925_m, 950_m, 975_m,
                                   1000_m};

    for (const auto &threshold : thresholds) {
        sims::Optics optics{};

        auto ordering = optics.compute_ordering(nodes, 0.0000001, MIN_PTS);
        auto clusters = optics.cluster(ordering);
        auto links = sims::data::create_link_vector(clusters, threshold);

        std::cout << "Distance threshold: " << threshold << std::endl;
        std::cout << "Links: " << links.size() << std::endl;

        if (links.size() <= 25000) {
            auto corr = sims::math::generate_correlation_matrix(links);
            auto std_deviation = std::pow(STANDARD_DEVIATION, 2);
            auto sigma = std_deviation * corr;

            auto chol_start = std::chrono::high_resolution_clock::now();
            auto autocorrelation_matrix = sims::cholesky::cholesky_v2(sigma);
            auto chol_time = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - chol_start);

            std::cout << "Cholesky time: " << format_duration(chol_time) << std::endl;
            std::cout << "\n" << std::endl;
        } else {
            std::cout << "Cholesky time: NaN\n\n" << std::endl;

        }
    }
}

#if 0 /* This block is for clustering experiments on eps */
int main(int argc, char *argv[]) {
    geo::Location upper{57.01266813458001, 10.994625734716218};
    auto lower = geo::square(upper, 10_km);
    auto nodes = sims::data::generate_nodes(NODES, upper, lower);

    std::vector<int> cluster_sizes{1000, 900, 800, 700, 600, 500, 450, 400, 350, 300, 250, 200, 150, 100, 50};
    for (const auto &k : cluster_sizes) {
        auto conv = true;
        auto eps = MIN_EPS;

        sims::Optics optics{};
        std::vector<sims::Optics::Cluster> clusters;

        std::cout << "Nodes: " << NODES << std::endl;
        while (conv) {
            optics = sims::Optics{};

            auto cluster_start = std::chrono::high_resolution_clock::now();
            auto ordering = optics.compute_ordering(nodes, eps, MIN_PTS);
            clusters = optics.cluster(ordering);
            auto cluster_time = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - cluster_start);

            if (clusters.size() <= k) {
                std::cout << "Clustering time: " << format_duration(cluster_time) << std::endl;
                std::cout << "Cluster goal: " << k << std::endl;
                std::cout << "Actual cluster count: " << clusters.size() << std::endl;
                std::cout << "eps: " << eps << std::endl;

                conv = false;
                eps = MIN_EPS;
            } else {
                eps += EPS_STEP;
            }
        }

        auto links = sims::data::create_link_vector(clusters);
        std::cout << "Links: " << links.size() << std::endl;
        if (links.size() <= 25000) {
            auto corr = sims::math::generate_correlation_matrix(links);
            auto std_deviation = std::pow(STANDARD_DEVIATION, 2);
            auto sigma = std_deviation * corr;

            auto chol_start = std::chrono::high_resolution_clock::now();
            auto autocorrelation_matrix = sims::cholesky::cholesky_v2(sigma);
            auto chol_time = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - chol_start);

            std::cout << "Cholesky time: " << format_duration(chol_time) << std::endl;
            std::cout << "\n" << std::endl;
        } else {
            std::cout << "Cholesky time: NaN\n\n" << std::endl;
        }

    }
}
#endif

#if 0
int main(int argc, char *argv[]) {
    geo::Location upper{57.01266813458001, 10.994625734716218};
    auto lower = geo::square(upper, 10_km);
    auto nodes = sims::data::generate_nodes(NODES, upper, lower);

    auto links = sims::data::create_link_vector(nodes);
    auto correlation_matrix = sims::math::generate_correlation_matrix(links);
    auto correlation_norm = sims::linalg::frobenius_norm(correlation_matrix);
    std::cout << "nodes: " << nodes.size() << std::endl;
    std::cout << "links: " << links.size() << std::endl;
    //std::cout << "norm: " << norm << std::endl;
    std::cout << std::endl;

    int eps_step_counter = 0;
    for (auto minpts = MIN_PTS; minpts <= MAX_PTS; minpts += PTS_STEPS) {
        for (auto eps = MIN_EPS; eps_step_counter < MAX_EPS_STEPS; eps += EPS_STEP, ++eps_step_counter) {
            std::cout << "eps: " << eps << std::endl;
            //std::cout << "minpts: " << minpts << std::endl;

            sims::Optics optics{};

            auto start = std::chrono::high_resolution_clock::now();

            auto ordering = optics.compute_ordering(nodes, eps, minpts);
            auto clusters = optics.cluster(ordering);
            //auto cluster_links = sims::data::create_link_vector(clusters);
            //auto clusters_correlation_matrix = sims::math::generate_correlation_matrix(clusters_links);
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                    std::chrono::high_resolution_clock::now() - start);

            std::vector<sims::Node> centroids{};
            centroids.reserve(nodes.size());
            for (auto &node : nodes) {
                for (auto &cluster : clusters) {
                    if (!cluster.contains(node)) {
                        continue;
                    }

                    centroids.emplace_back(node.get_id(), cluster.centroid());
                    break;
                }
            }

            auto centroid_links = sims::data::create_link_vector(centroids);
            auto centroid_correlations = sims::math::generate_correlation_matrix(centroid_links);
            auto corr_diff = correlation_matrix - centroid_correlations;
            auto norm = sims::linalg::frobenius_norm(corr_diff);
            auto centroid_norm = sims::linalg::frobenius_norm(centroid_correlations);

            std::cout << "clusters: " << clusters.size() << std::endl;
            //std::cout << "cluster_links: " << cluster_links.size() << std::endl;
            std::cout << "centroid_links: " << centroid_links.size() << std::endl;
            std::cout << "centroid_norm: " << centroid_norm << std::endl;
            std::cout << "norm_diff: " << centroid_norm - correlation_norm << std::endl;
            std::cout << "norm: " << norm << std::endl;
            std::cout << "duration: " << format_duration(duration) << std::endl;
            std::cout << std::endl;
        }
    }
}
#endif

#if 0

#define TRY try {

#define CATCH   } catch (const char * e) {           \
                    std::cout << e << std::endl;     \
                    return 0;                        \
                }


std::vector<sims::Optics::Cluster> cluster(std::vector<sims::Node> &nodes, double eps, int minpts) {
    sims::Optics optics{};

    auto ordering = optics.compute_ordering(nodes, eps, minpts);
    auto clusters = optics.cluster(ordering);

    return clusters;
}


void print_header() {
    std::cout
            << std::setw(12) << "Run"
            << std::setw(12) << "Clusters"
            << std::setw(12) << "Epsilon"
            << std::setw(12) << "Nodes"
            << std::setw(12) << "Links"
            << std::setw(12) << "Delta"
            << std::setw(12) << "Time"
            << std::endl;
}

void print_result(int graph_id, long cluster_count, double eps, long node_count, long link_count, double delta,
                  long duration) {
    double us = duration / 1000.0;
    long ms = static_cast<long>(std::round(us));
    std::cout
            << std::setprecision(3)
            << std::fixed
            << std::setw(8) << "" << std::setw(4) << std::setfill('0') << graph_id << std::setfill(' ')
            << std::setw(12) << cluster_count
            << std::setw(9) << eps << " km"
            << std::setw(12) << node_count
            << std::setw(12) << link_count
            << std::setw(12) << delta;

    if (ms > 100) {
        std::cout
                << std::setw(9) << ms << " ms";
    } else {
        std::cout
                << std::setw(9) << us << " ms";
    }
    std::cout
            << std::endl;
}
void request_graph(mpilib::HttpClient &httpclient, std::vector<sims::Optics::Cluster> clusters, int graph_id, double eps) {
    std::vector<json> serialized_clusters{};
    serialized_clusters.reserve(clusters.size());

    for (auto &cluster : clusters) {
        std::vector<int> node_ids{};
        node_ids.reserve(cluster.get_nodes().size());

        for (auto &node : cluster.get_nodes()) {
            node_ids.emplace_back(node.get_id());
        }

        json c = node_ids;
        serialized_clusters.emplace_back(c);
    }

    json json_clusters = {
            {"clusters", serialized_clusters},
            {"params",   {
                                 {"id", graph_id},
                                 {"minpts", MIN_PTS},
                                 {"eps", eps},
                                 {"count", clusters.size()}
                         }}
    };

    httpclient.post("/request-graph", json_clusters);
}

long clusterize_remaining(const std::vector<sims::Node> &nodes, std::vector<sims::Optics::Cluster> &clusters) {
    auto cnodes = nodes;
    unsigned long node_count = 0;
    for (auto &cluster : clusters) {
        node_count += cluster.size();

        for (auto &node : cluster.get_nodes()) {
            if (cluster.contains(node)) {
                cnodes.erase(std::remove(cnodes.begin(), cnodes.end(), node), cnodes.end());
            }
        }
    }
    node_count = nodes.size() - node_count + clusters.size();
    auto id = static_cast<uint32_t>(clusters.size());
    for (auto &node : cnodes) {
        id++;
        std::vector<sims::Node> single_node_cluster{node};
        sims::Optics::Cluster c{id, single_node_cluster};
        clusters.push_back(c);
    }

    assert(node_count == clusters.size());
    return node_count;
}

geo::Location create_square(geo::Location &upper, double size) {
    return move_location(move_location(upper, size, 180), size, 90);
}

int main(int argc, char *argv[]) {
    std::vector<std::future<void>> futures{};
    std::vector<std::thread> threads{};

    geo::Location upper{57.0134, 9.99008};
    geo::Location lower = create_square(upper, AREA);

    auto setup_start = std::chrono::steady_clock::now();
    /**/

    auto nodes = sims::data::generate_nodes(NODES, upper, lower);

    std::vector<sims::Optics::Cluster> og_clusters{};
    og_clusters.reserve(nodes.size());

    for (auto &node : nodes) {
        std::vector<sims::Node> cluster{node};
        og_clusters.emplace_back(node.get_id(), cluster);
    }
    auto og_links = sims::data::create_link_vector(og_clusters, LINK_THRESHOLD);
    std::cout << og_links.size() << std::endl;

    /* Compute distance part */
    std::vector<double> og_l_distance{};
    std::for_each(og_links.cbegin(), og_links.cend(), [&og_l_distance](sims::Optics::CLink link) {
        og_l_distance.emplace_back(sims::math::distance_pathloss(link.get_distance()));
    });

    /* Compute fading part */
    auto og_corr = sims::math::generate_correlation_matrix(og_links);
    std::cout << "pd: " << sims::cholesky::is_positive_definite(og_corr) << std::endl;

    if (!sims::cholesky::is_positive_definite(og_corr)) {
//        og_corr = iterative_spectral(og_corr, 25, 0.001);
        return 0;
    }

    std::cout << "pd: "  << sims::cholesky::is_positive_definite(og_corr) << std::endl;

//    auto setup_end1 = std::chrono::steady_clock::now();
//    auto setup_duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(setup_end1 - setup_start);


    auto std_deviation = std::pow(11.4, 2);
    auto og_sigma = std_deviation * og_corr;
    auto og_gaussian = sims::math::generate_gaussian_vector(0.0, 1.0, og_links.size());
    auto og_choleskied = sims::cholesky::cholesky_v2(og_sigma);
    assert(sims::linalg::compare_vectors(og_sigma, og_choleskied * sims::linalg::transpose(og_choleskied), 0.00001));
    auto og_l_fading = og_choleskied * og_gaussian;

    auto og_rssi = TX_DBM - (og_l_distance + og_l_fading);

    std::vector<double> og_pep{};
    og_pep.reserve(og_links.size());

    for (auto &rssi : og_rssi) {
        og_pep.push_back(sims::radiomodel::pep(rssi, 160));
    }

    /**/
    auto setup_end = std::chrono::steady_clock::now();
    auto setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(setup_end - setup_start);

    mpilib::HttpClient httpclient{"http://0.0.0.0:5000/vis"};
    json j_nodes = nodes;
    httpclient.post("/add-nodes", j_nodes);

    auto eps_steps = static_cast<unsigned long>(std::round(MAX_EPS / EPS_STEP));

    print_header();

    auto m = 0;

    print_result(0, 0, 0.0, nodes.size(), og_links.size(), 0.0, setup_duration.count());

    for (int i = 0; i < eps_steps; ++i) {
        m++;
        auto eps = MIN_EPS + i * EPS_STEP;

        auto start = std::chrono::steady_clock::now();
        /**/

        auto clusters = cluster(nodes, eps, MIN_PTS);

        auto cluster_count = clusters.size();
        auto node_count = clusterize_remaining(nodes, clusters);
        auto links = sims::data::create_link_vector(clusters, LINK_THRESHOLD);

        /* Compute distance part */
        std::vector<double> l_distance{};
        std::for_each(links.cbegin(), links.cend(), [&l_distance](sims::Optics::CLink link) {
            l_distance.emplace_back(sims::math::distance_pathloss(link.get_distance()));
        });

        /* Compute fading part */
        auto corr = sims::math::generate_correlation_matrix(links);
        auto sigma = std_deviation * corr;
        auto gaussian = sims::math::generate_gaussian_vector(0.0, 1.0, links.size());
        auto l_fading = sims::cholesky::cholesky(sigma) * gaussian;

        std::vector<double> rssi{TX_DBM - (l_distance + l_fading)};

        assert(rssi.size() == links.size());
        for (int j = 0; j < rssi.size(); ++j) {
            links[j].set_rssi(rssi[j]);
        }

        for (auto &og_link : og_links) {
            auto first = og_link.get_clusters().first.get_nodes().front();
            auto second = og_link.get_clusters().second.get_nodes().front();
            auto k = 0;
            for (auto &link : links) {
                k++;
                if (!link.contains(first, second)) {
                    continue;
                }
                link.set_pep(sims::radiomodel::pep(link.get_rssi(), 160));
                std::cout
                        << link.get_rssi()
                        << " : "
                        << link.get_pep()
                        << std::endl;
            }
        }

        /**/
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        print_result(m, cluster_count, eps, node_count, links.size(), 0.0, duration.count());

        if (clusters.empty()) {
            continue;
        }

        std::future<void> future = std::async(std::launch::async, request_graph,
                                              std::ref(httpclient), clusters, m, eps);
        futures.push_back(std::move(future));
    }

    for (auto &future : futures) {
        future.wait();
    }
}



#endif
