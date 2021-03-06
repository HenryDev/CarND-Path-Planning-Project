#include <fstream>
#include <cmath>
#include <uWS/uWS.h>
#include <thread>
#include "Eigen-3.3/Eigen/Core"
#include "json.hpp"
#include "spline.h"

using namespace std;

using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(const string &s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of('[');
    auto b2 = s.find_first_of('}');
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double>
getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y) {
    int prev_wp = -1;

    while (s > maps_s[prev_wp + 1] && (prev_wp < (int) (maps_s.size() - 1))) {
        prev_wp++;
    }

    auto wp2 = static_cast<int>((prev_wp + 1) % maps_x.size());

    double heading = atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s - maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
    double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

    double perp_heading = heading - pi() / 2;

    double x = seg_x + d * cos(perp_heading);
    double y = seg_y + d * sin(perp_heading);

    return {x, y};

}

int main() {
    uWS::Hub h;

    // Load up map values for waypoint's x,y,s and d normalized normal vectors
    vector<double> map_waypoints_x;
    vector<double> map_waypoints_y;
    vector<double> map_waypoints_s;
    vector<double> map_waypoints_dx;
    vector<double> map_waypoints_dy;

    // Waypoint map to read from
    string map_file_ = "../data/highway_map.csv";
    // The max s value before wrapping around the track back to 0
    double max_s = 6945.554;

    ifstream in_map_(map_file_.c_str(), ifstream::in);

    string line;
    while (getline(in_map_, line)) {
        istringstream iss(line);
        double x;
        double y;
        double s;
        double d_x;
        double d_y;
        iss >> x;
        iss >> y;
        iss >> s;
        iss >> d_x;
        iss >> d_y;
        map_waypoints_x.push_back(x);
        map_waypoints_y.push_back(y);
        map_waypoints_s.push_back(s);
        map_waypoints_dx.push_back(d_x);
        map_waypoints_dy.push_back(d_y);
    }
    int lane = 1; //start in lane 1 the middle lane
    double reference_velocity = 0; //starting speed in mph
    h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s, &map_waypoints_dx, &map_waypoints_dy,
                        &lane, &reference_velocity](
            uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
            uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        //auto sdata = string(data).substr(0, length);
        //cout << sdata << endl;
        if (length > 2 && data[0] == '4' && data[1] == '2') {

            auto basicString = hasData(data);

            if (!basicString.empty()) {
                auto j = json::parse(basicString);

                string event = j[0].get<string>();

                if (event == "telemetry") {
                    // j[1] is the data JSON object

                    // Main car's localization Data
                    double car_x = j[1]["x"];
                    double car_y = j[1]["y"];
                    double car_s = j[1]["s"];
                    double car_d = j[1]["d"];
                    double car_yaw = j[1]["yaw"];
                    double car_speed = j[1]["speed"];

                    // Previous path data given to the Planner
                    auto previous_path_x = j[1]["previous_path_x"];
                    auto previous_path_y = j[1]["previous_path_y"];
                    // Previous path's end s and d values
                    double end_path_s = j[1]["end_path_s"];
                    double end_path_d = j[1]["end_path_d"];

                    // Sensor Fusion Data, a list of all other cars on the same side of the road.
                    vector<vector<double>> sensor_fusion = j[1]["sensor_fusion"];
                    //previous list of points
                    int previous_size = previous_path_x.size();

                    if (previous_size > 0) {
                        car_s = end_path_s;
                    }
                    bool too_close = false;
                    bool middle_lane_occupied = false;
                    bool left_lane_occupied = false;
                    bool right_lane_occupied = false;
                    double breaking_magnitude = 0.5;
                    for (auto &sensor_data : sensor_fusion) {
                        //the s value of that car
                        double their_s = sensor_data[5];
                        //car's d position in frenet coordinates
                        double their_d = sensor_data[6];
                        double vx = sensor_data[3];
                        double vy = sensor_data[4];
                        double their_speed = sqrt(vx * vx + vy * vy);
                        //if using previous path points, project the s values outward in time
                        their_s += previous_size * 0.02 * their_speed;

                        if ((their_s > car_s && their_s - car_s < 35) || (their_s < car_s && car_s - their_s < 10)) {
                            if (their_d <= 4) {
                                left_lane_occupied = true;
                            }
                            if (their_d > 4 && their_d <= 8) {
                                middle_lane_occupied = true;
                            }
                            if (their_d > 8) {
                                right_lane_occupied = true;
                            }
                        }

                        //if the car is in my lane
                        if (their_d < 2 + 4 * lane + 2 && their_d > 2 + 4 * lane - 2) {
                            //if they're in front of us and the gap is less than 30m
                            if (their_s > car_s && their_s - car_s < 20) {
                                too_close = true;
                                if (left_lane_occupied && middle_lane_occupied && right_lane_occupied) {
                                    reference_velocity -= breaking_magnitude;
                                    continue;
                                }
                                if (lane == 0 || lane == 2) {
                                    if (!middle_lane_occupied) {
                                        lane = 1;
                                    } else {
                                        reference_velocity -= breaking_magnitude;
                                    }
                                } else if (lane == 1) {
                                    if (!left_lane_occupied) {
                                        lane = 0;
                                    } else if (!right_lane_occupied) {
                                        lane = 2;
                                    } else {
                                        reference_velocity -= breaking_magnitude;
                                    }
                                }
                            }
                        }
                    }
                    if (too_close) {
                        //slow down gradually
                        reference_velocity -= breaking_magnitude;
                    } else if (reference_velocity < 49) {
                        //speed up gradually
                        reference_velocity += 1;
                    }

                    // waypoints
                    vector<double> ptsx;
                    vector<double> ptsy;

                    //where the car is or where the previous path ends
                    double reference_x = car_x;
                    double reference_y = car_y;
                    double reference_yaw = deg2rad(car_yaw);
                    //if the previous path is almost empty, then use the car state as a starting reference
                    if (previous_size < 2) {
                        //use 2 points that are tangent to where the car is pointing
                        double previous_car_x = car_x - cos(car_yaw);
                        double previous_car_y = car_y - sin(car_yaw);
                        //generate 2 points based on the car's angle
                        ptsx.push_back(previous_car_x);
                        ptsx.push_back(car_x);
                        ptsy.push_back(previous_car_y);
                        ptsy.push_back(car_y);
                    } else {
                        //use the last 2 points in the previous path to calculate where the car is heading
                        reference_x = previous_path_x[previous_size - 1];
                        reference_y = previous_path_y[previous_size - 1];
                        double reference_previous_x = previous_path_x[previous_size - 2];
                        double reference_previous_y = previous_path_y[previous_size - 2];
                        reference_yaw = atan2(reference_y - reference_previous_y, reference_x - reference_previous_x);
                        ptsx.push_back(reference_previous_x);
                        ptsx.push_back(reference_x);
                        ptsy.push_back(reference_previous_y);
                        ptsy.push_back(reference_y);
                    }

                    // use spline to smooth lane shifts
                    int d = 2 + 4 * lane;
                    double s = car_s + 30;
                    double s1 = car_s + 60;
                    double s2 = car_s + 90;
                    //make 3 points 30m apart for the location of the car
                    const vector<double> &next_wp0 = getXY(s, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    const vector<double> &next_wp1 = getXY(s1, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    const vector<double> &next_wp2 = getXY(s2, d, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    ptsx.push_back(next_wp0[0]);
                    ptsx.push_back(next_wp1[0]);
                    ptsx.push_back(next_wp2[0]);
                    ptsy.push_back(next_wp0[1]);
                    ptsy.push_back(next_wp1[1]);
                    ptsy.push_back(next_wp2[1]);

                    //shift the perspective so that the points' heading is the same as the car's heading
                    for (int i = 0; i < ptsx.size(); i++) {
                        double shift_x = ptsx[i] - reference_x;
                        double shift_y = ptsy[i] - reference_y;
                        ptsx[i] = shift_x * cos(-reference_yaw) - shift_y * sin(-reference_yaw);
                        ptsy[i] = shift_x * sin(-reference_yaw) + shift_y * cos(-reference_yaw);
                    }

                    tk::spline spline;
                    //add the 5 anchor points to spline
                    spline.set_points(ptsx, ptsy);

                    //the actual x and y points for use in the planner
                    vector<double> next_x_vals;
                    vector<double> next_y_vals;
                    //if there's any points leftover from the previous path, add them to the path planner
                    for (int i = 0; i < previous_path_x.size(); i++) {
                        next_x_vals.push_back(previous_path_x[i]);
                        next_y_vals.push_back(previous_path_y[i]);
                    }
                    //horizon is 30m
                    double target_x = 30.0;
                    double target_y = spline(target_x);
                    //from car to target
                    double target_distance = sqrt(target_x * target_x + target_y * target_y);
                    //where we start at the begining of this calculation
                    double x_add_on = 0;
                    for (int i = 1; i < 50 - previous_path_x.size(); i++) {
                        // N * 0.02 * velocity = d
                        double n = target_distance / (0.02 * reference_velocity / 2.24); //mph to m/s
                        double x_point = x_add_on + target_x / n;
                        double y_point = spline(x_point);
                        x_add_on = x_point;
                        double x_reference = x_point;
                        double y_reference = y_point;
                        //shift and rotate back to map coordinates
                        x_point = x_reference * cos(reference_yaw) - y_reference * sin(reference_yaw);
                        y_point = x_reference * sin(reference_yaw) + y_reference * cos(reference_yaw);
                        x_point += reference_x;
                        y_point += reference_y;
                        next_x_vals.push_back(x_point);
                        next_y_vals.push_back(y_point);
                    }

                    json msgJson;
                    msgJson["next_x"] = next_x_vals;
                    msgJson["next_y"] = next_y_vals;

                    auto msg = "42[\"control\"," + msgJson.dump() + "]";

                    //this_thread::sleep_for(chrono::milliseconds(1000));
                    ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                }
            } else {
                // Manual driving
                std::string msg = "42[\"manual\",{}]";
                ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
        }
    });

    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
}