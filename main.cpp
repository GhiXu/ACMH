#include "main.h"
#include "ACMH.h"

void GenerateSampleList(const std::string &dense_folder, std::vector<Problem> &problems)
{
    std::string cluster_list_path = dense_folder + std::string("/pair.txt");

    problems.clear();

    std::ifstream file(cluster_list_path);

    int num_images;
    file >> num_images;

    for (int i = 0; i < num_images; ++i) {
        Problem problem;
        problem.src_image_ids.clear();
        file >> problem.ref_image_id;

        int num_src_images;
        file >> num_src_images;
        for (int j = 0; j < num_src_images; ++j) {
            int id;
            float score;
            file >> id >> score;
            if (score <= 0.0f) {
                continue;
            }
            problem.src_image_ids.push_back(id);
        }
        problems.push_back(problem);
    }
}

void ProcessProblem(const std::string &dense_folder, const Problem &problem, bool geom_consistency)
{
    std::cout << "Processing image " << std::setw(8) << std::setfill('0') << problem.ref_image_id << "..." << std::endl;
    cudaSetDevice(1);
    std::stringstream result_path;
    result_path << dense_folder << "/ACMH" << "/2333_" << std::setw(8) << std::setfill('0') << problem.ref_image_id;
    std::string result_folder = result_path.str();
    mkdir(result_folder.c_str(), 0777);

    ACMH acmh;
    if (geom_consistency) {
        acmh.SetGeomConsistencyParams();
    }
    acmh.InuputInitialization(dense_folder, problem);

    acmh.CudaSpaceInitialization(dense_folder, problem);
    acmh.RunPatchMatch();

    const int width = acmh.GetReferenceImageWidth();
    const int height = acmh.GetReferenceImageHeight();

    cv::Mat_<float> depths = cv::Mat::zeros(height, width, CV_32FC1);
    cv::Mat_<cv::Vec3f> normals = cv::Mat::zeros(height, width, CV_32FC3);
    cv::Mat_<float> costs = cv::Mat::zeros(height, width, CV_32FC1);

    for (int col = 0; col < width; ++col) {
        for (int row = 0; row < height; ++row) {
            int center = row * width + col;
            float4 plane_hypothesis = acmh.GetPlaneHypothesis(center);
            depths(row, col) = plane_hypothesis.w;
            normals(row, col) = cv::Vec3f(plane_hypothesis.x, plane_hypothesis.y, plane_hypothesis.z);
            costs(row, col) = acmh.GetCost(center);
        }
    }

    std::string suffix = "/depths.dmb";
    if (geom_consistency) {
        suffix = "/depths_geom.dmb";
    }
    std::string depth_path = result_folder + suffix;
    std::string normal_path = result_folder + "/normals.dmb";
    std::string cost_path = result_folder + "/costs.dmb";
    writeDepthDmb(depth_path, depths);
    writeNormalDmb(normal_path, normals);
    writeDepthDmb(cost_path, costs);
    std::cout << "Processing image " << std::setw(8) << std::setfill('0') << problem.ref_image_id << " done!" << std::endl;
}

void RunFusion(std::string &dense_folder, const std::vector<Problem> &problems, bool geom_consistency)
{
    size_t num_images = problems.size();
    std::string image_folder = dense_folder + std::string("/images");
    std::string cam_folder = dense_folder + std::string("/cams");

    std::vector<cv::Mat> images;
    std::vector<Camera> cameras;
    std::vector<cv::Mat_<float>> depths;
    std::vector<cv::Mat_<cv::Vec3f>> normals;
    std::vector<cv::Mat> masks;
    images.clear();
    cameras.clear();
    depths.clear();
    normals.clear();
    masks.clear();

    for (size_t i = 0; i < num_images; ++i) {
        std::cout << "Reading image " << std::setw(8) << std::setfill('0') << i << "..." << std::endl;
        std::stringstream image_path;
        image_path << image_folder << "/" << std::setw(8) << std::setfill('0') << problems[i].ref_image_id << ".jpg";
        cv::Mat_<cv::Vec3b> image = cv::imread (image_path.str(), cv::IMREAD_COLOR);
        std::stringstream cam_path;
        cam_path << cam_folder << "/" << std::setw(8) << std::setfill('0') << problems[i].ref_image_id << "_cam.txt";
        Camera camera = ReadCamera(cam_path.str());

        std::stringstream result_path;
        result_path << dense_folder << "/ACMH" << "/2333_" << std::setw(8) << std::setfill('0') << problems[i].ref_image_id;
        std::string result_folder = result_path.str();
        std::string suffix = "/depths.dmb";
        if (geom_consistency) {
            suffix = "/depths_geom.dmb";
        }
        std::string depth_path = result_folder + suffix;
        std::string normal_path = result_folder + "/normals.dmb";
        cv::Mat_<float> depth;
        cv::Mat_<cv::Vec3f> normal;
        readDepthDmb(depth_path, depth);
        readNormalDmb(normal_path, normal);

        cv::Mat_<cv::Vec3b> scaled_image;
        RescaleImageAndCamera(image, scaled_image, depth, camera);
        images.push_back(scaled_image);
        cameras.push_back(camera);
        depths.push_back(depth);
        normals.push_back(normal);
        cv::Mat mask = cv::Mat::zeros(depth.rows, depth.cols, CV_8UC1);
        masks.push_back(mask);
    }

    std::vector<PointList> PointCloud;
    PointCloud.clear();

    for (size_t i = 0; i < num_images; ++i) {
        std::cout << "Fusing image " << std::setw(8) << std::setfill('0') << i << "..." << std::endl;
        const int cols = depths[i].cols;
        const int rows = depths[i].rows;
        int num_ngb = problems[i].src_image_ids.size();
        std::vector<int2> used_list(num_ngb, make_int2(-1, -1));
        for (int r =0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (masks[i].at<uchar>(r, c) == 1)
                    continue;
                float ref_depth = depths[i].at<float>(r, c);
                cv::Vec3f ref_normal = normals[i].at<cv::Vec3f>(r, c);

                if (ref_depth <= 0.0)
                    continue;

                float3 PointX = Get3DPointonWorld(c, r, ref_depth, cameras[i]);
                float3 consistent_Point = PointX;
                cv::Vec3f consistent_normal = ref_normal;
                float consistent_Color[3] = {(float)images[i].at<cv::Vec3b>(r, c)[0], (float)images[i].at<cv::Vec3b>(r, c)[1], (float)images[i].at<cv::Vec3b>(r, c)[2]};
                int num_consistent = 0;

                for (int j = 0; j < num_ngb; ++j) {
                    int src_id = problems[i].src_image_ids[j];
                    const int src_cols = depths[src_id].cols;
                    const int src_rows = depths[src_id].rows;
                    float2 point;
                    float proj_depth;
                    ProjectonCamera(PointX, cameras[src_id], point, proj_depth);
                    int src_r = int(point.y + 0.5f);
                    int src_c = int(point.x + 0.5f);
                    if (src_c >= 0 && src_c < src_cols && src_r >= 0 && src_r < src_rows) {
                        if (masks[src_id].at<uchar>(src_r, src_c) == 1)
                            continue;

                        float src_depth = depths[src_id].at<float>(src_r, src_c);
                        cv::Vec3f src_normal = normals[src_id].at<cv::Vec3f>(src_r, src_c);
                        if (src_depth <= 0.0)
                            continue;

                        float3 tmp_X = Get3DPointonWorld(src_c, src_r, src_depth, cameras[src_id]);
                        float2 tmp_pt;
                        ProjectonCamera(tmp_X, cameras[i], tmp_pt, proj_depth);
                        float reproj_error = sqrt(pow(c - tmp_pt.x, 2) + pow(r - tmp_pt.y, 2));
                        float relative_depth_diff = fabs(proj_depth - ref_depth) / ref_depth;
                        float angle = GetAngle(ref_normal, src_normal);

                        if (reproj_error < 2.0f && relative_depth_diff < 0.01f && angle < 0.174533f) {
                            consistent_Point.x += tmp_X.x;
                            consistent_Point.y += tmp_X.y;
                            consistent_Point.z += tmp_X.z;
                            consistent_normal = consistent_normal + src_normal;
                            consistent_Color[0] += images[src_id].at<cv::Vec3b>(src_r, src_c)[0];
                            consistent_Color[1] += images[src_id].at<cv::Vec3b>(src_r, src_c)[1];
                            consistent_Color[2] += images[src_id].at<cv::Vec3b>(src_r, src_c)[2];

                            used_list[j].x = src_c;
                            used_list[j].y = src_r;
                            num_consistent++;
                        }
                    }
                }

                if (num_consistent >= 2) {
                    consistent_Point.x /= (num_consistent + 1.0f);
                    consistent_Point.y /= (num_consistent + 1.0f);
                    consistent_Point.z /= (num_consistent + 1.0f);
                    consistent_normal /= (num_consistent + 1.0f);
                    consistent_Color[0] /= (num_consistent + 1.0f);
                    consistent_Color[1] /= (num_consistent + 1.0f);
                    consistent_Color[2] /= (num_consistent + 1.0f);

                    PointList point3D;
                    point3D.coord = consistent_Point;
                    point3D.normal = make_float3(consistent_normal[0], consistent_normal[1], consistent_normal[2]);
                    point3D.color = make_float3(consistent_Color[0], consistent_Color[1], consistent_Color[2]);
                    PointCloud.push_back(point3D);

                    for (int j = 0; j < num_ngb; ++j) {
                        if (used_list[j].x == -1)
                            continue;
                        masks[problems[i].src_image_ids[j]].at<uchar>(used_list[j].y, used_list[j].x) = 1;
                    }
                }
            }
        }
    }

    std::string ply_path = dense_folder + "/ACMH/ACMH_model.ply";
    StoreColorPlyFileBinaryPointCloud (ply_path, PointCloud);
}

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "USAGE: ACMH dense_folder" << std::endl;
        return -1;
    }

    std::string dense_folder = argv[1];
    std::vector<Problem> problems;
    GenerateSampleList(dense_folder, problems);

    std::string output_folder = dense_folder + std::string("/ACMH");
    mkdir(output_folder.c_str(), 0777);

    size_t num_images = problems.size();
    std::cout << "There are " << num_images << " problems needed to be processed!" << std::endl;

    bool geom_consistency = false;
    for (size_t i = 0; i < num_images; ++i) {
        ProcessProblem(dense_folder, problems[i], geom_consistency);
    }
    geom_consistency = true;
    for (size_t i = 0; i < num_images; ++i) {
        ProcessProblem(dense_folder, problems[i], geom_consistency);
    }

    RunFusion(dense_folder, problems, geom_consistency);

    return 0;
}
