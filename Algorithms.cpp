#include <iostream>
#include <vector>
#include <utility>
#include <queue>
#include <fstream>
#include <set>
#include <time.h>
#include <iomanip>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>

using namespace std;

const double earth_radius_km = 6371.0;

double deg2rad(double deg)
{
	return (deg * M_PI / 180.0);
}

double haversine_distance(double latitude1, double longitude1, double latitude2,
	double longitude2)
{
	double lat1 = deg2rad(latitude1);
	double lon1 = deg2rad(longitude1);
	double lat2 = deg2rad(latitude2);
	double lon2 = deg2rad(longitude2);

	double d_lat = abs(lat1 - lat2);
	double d_lon = abs(lon1 - lon2);

	double a = pow(sin(d_lat / 2), 2) + cos(lat1) * cos(lat2) * pow(sin(d_lon / 2), 2);

	//double d_sigma = 2 * atan2(sqrt(a), sqrt(1 - a));
	double d_sigma = 2 * asin(sqrt(a));

	return earth_radius_km * d_sigma;
}

struct intersection {
	double lon, lat;

	intersection() {
		lon = 0;
		lat = 0;
	}

	intersection(double x, double y) {
		lon = x;
		lat = y;
	}

};

pair<vector<int>,vector<pair<int,int>>> dijkstra(vector<vector<pair<double, int>>> &g, vector<intersection> &v, int s, int e) {
	
	priority_queue<pair<double, int>> pq;
	vector<double> dist(v.size(), DBL_MAX);
	vector<int> prev(v.size(), -1);
	vector<pair<int, int>> checked;

	dist[s] = 0;
	pq.push({ -dist[s], s });

	while (pq.top().second != e) {

		pair<double, int> cur = pq.top();
		pq.pop();
		int cur_vertex = cur.second;
		double cur_dist = -cur.first;
		if (cur_dist > dist[cur_vertex]) {
			continue;
		}
		
		for (int i = 0; i < g[cur_vertex].size(); i++) {
			int to = g[cur_vertex][i].second;
			double to_dist = g[cur_vertex][i].first;

			checked.push_back(make_pair(cur_vertex, to));
				
			if (cur_dist + to_dist < dist[to]) {
				dist[to] = cur_dist + to_dist;
				prev[to] = cur_vertex;
				pq.push(make_pair(-dist[to], to));
			}
		}

		if (pq.empty()) {
			return { {}, checked };
		}

	}
	vector<int> res;
	int id = e;
	while (id != -1) {
		res.push_back(id);
		id = prev[id];
	}
	reverse(res.begin(), res.end());
	cout << "dijkstra res " << dist[e] << endl;
	return { res, checked };

}

pair<vector<int>, vector<pair<int, int>>> biderectional_dijkstra(vector<vector<pair<double, int>>>& g, vector<intersection>& v, int s, int e) {

	double inf = 1e10;

	priority_queue<pair<double, int>> pq_forward;
	priority_queue<pair<double, int>> pq_backward;
	vector<double> dist_forward(v.size(), inf);
	vector<double> dist_backward(v.size(), inf);
	vector<int> prev_forward(v.size(), -1);
	vector<int> prev_backward(v.size(), -1);
	vector<pair<int, int>> checked;
	vector<bool> poped_forward(v.size());
	vector<bool> poped_backward(v.size());

	double min_dist = inf;

	dist_forward[s] = 0;
	dist_backward[e] = 0;
	pq_forward.push({ -dist_forward[s], s });
	pq_backward.push({ -dist_backward[e], e });

	int iter = 0;

	int connection_vertex = -2;

	while (true) {

		
		pair<double, int> cur_forward = pq_forward.top();
		

		int cur_vertex_forward = cur_forward.second;
		double cur_dist_forward = -cur_forward.first;
		while (cur_dist_forward > dist_forward[cur_vertex_forward]) {
			pq_forward.pop();
			cur_forward = pq_forward.top();
			cur_vertex_forward = cur_forward.second;
			cur_dist_forward = -cur_forward.first;
		}

		pair<double, int> cur_backward = pq_backward.top();

		int cur_vertex_backward = cur_backward.second;
		double cur_dist_backward = -cur_backward.first;
		while (cur_dist_backward > dist_backward[cur_vertex_backward]) {
			pq_backward.pop();
			cur_backward = pq_backward.top();
			cur_vertex_backward = cur_backward.second;
			cur_dist_backward = -cur_backward.first;
			
		}

		if (dist_forward[cur_vertex_forward] + dist_backward[cur_vertex_backward] >= min_dist) {
			break;
		}

		if (cur_dist_forward < cur_dist_backward) {
			
			for (int i = 0; i < g[cur_vertex_forward].size(); i++) {
				int to = g[cur_vertex_forward][i].second;
				double to_dist = g[cur_vertex_forward][i].first;

				checked.push_back(make_pair(cur_vertex_forward, to));

				if (cur_dist_forward + to_dist < dist_forward[to]) {
					dist_forward[to] = cur_dist_forward + to_dist;
					prev_forward[to] = cur_vertex_forward;
					pq_forward.push(make_pair(-dist_forward[to], to));
				}

				if (poped_backward[to] && dist_forward[cur_vertex_forward] + to_dist + dist_backward[to] < min_dist) {
					min_dist = dist_forward[cur_vertex_forward] + to_dist + dist_backward[to];
					connection_vertex = to;
				}


			}



			poped_forward[cur_vertex_forward] = true;
			pq_forward.pop();
		}
		else {

			for (int i = 0; i < g[cur_vertex_backward].size(); i++) {
				int to = g[cur_vertex_backward][i].second;
				double to_dist = g[cur_vertex_backward][i].first;

				checked.push_back(make_pair(cur_vertex_backward, to));

				if (cur_dist_backward + to_dist < dist_backward[to]) {
					dist_backward[to] = cur_dist_backward + to_dist;
					prev_backward[to] = cur_vertex_backward;
					pq_backward.push(make_pair(-dist_backward[to], to));
				}

				if (poped_forward[to] && dist_backward[cur_vertex_backward] + to_dist + dist_forward[to] < min_dist) {
					min_dist = dist_backward[cur_vertex_backward] + to_dist + dist_forward[to];
					connection_vertex = to;
				}


			}
			poped_backward[cur_vertex_backward] = true;
			pq_backward.pop();
		}

		if (pq_forward.empty() || pq_backward.empty()) {
			return { {}, checked };
		}

	}

	

	if (poped_backward[pq_forward.top().second]) {
		connection_vertex = pq_forward.top().second;
	}
	else if(poped_forward[pq_backward.top().second]) {
		connection_vertex = pq_backward.top().second;
	}

	vector<int> res;
	int id = connection_vertex;
	while (id != -1) {
		res.push_back(id);
		id = prev_forward[id];
	}
	reverse(res.begin(), res.end());

	id = prev_backward[connection_vertex];
	while (id != -1) {
		res.push_back(id);
		id = prev_backward[id];
	}
	cout << "bidir res " << dist_forward[connection_vertex] + dist_backward[connection_vertex] << endl;
	return { res, checked };

}

double h(intersection x, intersection y) {
	return haversine_distance(x.lat, x.lon, y.lat, y.lon);
}

pair<vector<int>, vector<pair<int, int>>> astar(vector<vector<pair<double, int>>> &g, vector<intersection> &v, int s, int e) {

	priority_queue<pair<double, int>> pq;
	vector<double> g_score(v.size(), DBL_MAX);
	vector<double> f_score(v.size(), DBL_MAX);
	vector<int> prev(v.size(), -1);
	vector<pair<int, int>> checked;

	g_score[s] = 0;
	f_score[s] = h(v[e], v[s]);
	pq.push({ -h(v[e], v[s]), s });


	while (pq.top().second != e) {
		pair<double, int> cur = pq.top();
		pq.pop();
		int cur_vertex = cur.second;
		double cur_dist = -cur.first;
		if (cur_dist > f_score[cur_vertex]) {
			continue;
		}

		for (int i = 0; i < g[cur_vertex].size(); i++) {
			int to = g[cur_vertex][i].second;
			double to_dist = g[cur_vertex][i].first;

			checked.push_back(make_pair(cur_vertex, to));

			if (g_score[cur_vertex] + to_dist < g_score[to]) {
				g_score[to] = g_score[cur_vertex] + to_dist;
				f_score[to] = g_score[to] + h(v[e], v[to]);
				prev[to] = cur_vertex;
				pq.push(make_pair(-f_score[to], to));
			}
			
		}

		if (pq.empty()) {
			return { {}, checked };
		}

	}

	vector<int> res;
	int id = e;
	cout << "astar res " << g_score[e] << endl;
	while (id != -1) {
		res.push_back(id);
		id = prev[id];
	}
	reverse(res.begin(), res.end());
	return { res, checked };

}

void output_res(pair<vector<int>, vector<pair<int, int>>> ans, string filename) {
	vector<int> res = ans.first;
	vector<pair<int, int>> checked = ans.second;

	ofstream fout("C:\\Users\\aaa\\graphs\\"+filename);
	fout << res.size() << endl;
	for (int i = 1; i < res.size(); i++) {
		fout << res[i - 1] << " " << res[i] << endl;
	}
	fout << checked.size() << endl;
	for (int i = 0; i < checked.size(); i++) {
		fout << checked[i].first << " " << checked[i].second << endl;
	}
	fout.close();
}

int main() {

	string input_name = "input.txt";

	ifstream fin("C:\\Users\\aaa\\graphs\\"+ input_name);


	int n, m, s, e;
	fin >> n >> m >> s >> e;

	vector<vector<pair<double, int>>> g(n);
	vector<intersection> v(n);

	for (int i = 0; i < m; i++) {
		int u, v;
		double d;
		fin >> u >> v >> d;

		g[u].push_back(make_pair(d, v));
		g[v].push_back(make_pair(d, u));

	}

	for (int i = 0; i < n; i++) {
		double lon, lat;
		fin >> lon >> lat;
		v[i] = intersection(lon, lat);
	}

	fin.close();

	auto start_time = high_resolution_clock::now();
	auto ans = dijkstra(g, v, s, e);
	//output_res(ans, "output_dijkstra.txt");
	auto end_time = high_resolution_clock::now();
	auto duration = duration_cast<nanoseconds>(end_time - start_time);
	//cout << fixed << setprecision(7) << "dijkstra " << (double)duration.count() / 1e9 << endl << endl;
		


	start_time = high_resolution_clock::now();
	ans = biderectional_dijkstra(g, v, s, e);
	output_res(ans, "output_biderectional_dijkstra.txt");
	end_time = high_resolution_clock::now();
	duration = duration_cast<nanoseconds>(end_time - start_time);
	cout << fixed << setprecision(7) << "biderectional " << (double)duration.count() / 1e9 << endl << endl;

	start_time = high_resolution_clock::now();
	ans = astar(g, v, s, e);
	output_res(ans, "output_astar.txt");
	end_time = high_resolution_clock::now();
	duration = duration_cast<nanoseconds>(end_time - start_time);
	cout << fixed << setprecision(7) << "astar " << (double)duration.count() / 1e9 << endl << endl;

}