/*
   Taken from:
   The Computer Language Benchmarks Game
   https://salsa.debian.org/benchmarksgame-team/benchmarksgame/

   An implementation pretty much from scratch, with inspiration from the Rust
   version, which used the idea of saving some of the ingredients of the
   compution in an array instead of recomputing them.
   
   contributed by cvergu
   slightly modified by bmmeijers
*/

#define _USE_MATH_DEFINES // https://docs.microsoft.com/en-us/cpp/c-runtime-library/math-constants?view=msvc-160
#include <cmath>
#include <iostream>
#include <time.h>
#include <fstream>


// these values are constant and not allowed to be changed
const double SOLAR_MASS = 4 * M_PI * M_PI;
const double DAYS_PER_YEAR = 365.24;
const unsigned int BODIES_COUNT = 5;


class vector3d {
public:
    double x, y, z;

    double norm() const noexcept {
        return x * x + y * y + z * z;
    }

    double magnitude(double dt) const noexcept {
        double sum = norm();
        return dt / (sum * sqrt(sum));
    }
};

vector3d operator+(vector3d v1, vector3d v2) {
    return vector3d{
            v1.x + v2.x, v1.y + v2.y, v1.z + v2.z
    };
}

vector3d operator-(vector3d v1, vector3d v2) {
    return vector3d{
            v1.x - v2.x, v1.y - v2.y, v1.z - v2.z
    };
}

vector3d &operator+=(vector3d &v1, vector3d v2) {
    v1.x += v2.x;
    v1.y += v2.y;
    v1.z += v2.z;

    return v1;
}

vector3d &operator-=(vector3d &v1, vector3d v2) {
    v1.x -= v2.x;
    v1.y -= v2.y;
    v1.z -= v2.z;

    return v1;
}

vector3d &operator*=(vector3d &v, double mag) {
    v.x *= mag;
    v.y *= mag;
    v.z *= mag;

    return v;
}

vector3d operator*(vector3d v, double mag) {
    return vector3d{
            v.x * mag, v.y * mag, v.z * mag
    };
}

vector3d operator/(vector3d v, double mag) {
    return vector3d{
            v.x / mag, v.y / mag, v.z / mag
    };
}


class body {
public:
    std::string name;
    vector3d position;
    vector3d velocity;
    double mass;
};


void advance(body state[BODIES_COUNT], double dt) {
    /*
     * We precompute the quantity (r_i - r_j)
     */
    // 2D array (to hold: BODIES_COUNT x BODIES_COUNT elements)
    vector3d rij[BODIES_COUNT][BODIES_COUNT];

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            rij[i][j] = state[i].position - state[j].position;
        }
    }

    double magnitudes[BODIES_COUNT][BODIES_COUNT];

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            magnitudes[i][j] = rij[i][j].magnitude(dt);
        }
    }

    /*
     * Compute the new speed using v_i = a_i dt, where
     * a_i = \sum_{j \neq i} m_j (r_i - r_j)/|r_i - r_j|^3
     */
    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            vector3d dist = rij[i][j];
            double mag = magnitudes[i][j];
            state[i].velocity -= dist * (state[j].mass * mag);
            state[j].velocity += dist * (state[i].mass * mag);
        }
    }

    /*
     * Compute the new positions
     */
    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        state[i].position += state[i].velocity * dt;
    }
}

void offset_momentum(body state[BODIES_COUNT]) {
    vector3d &sun_velocity = state[0].velocity;

    for (unsigned int i = 1; i < BODIES_COUNT; ++i) {
        sun_velocity -= state[i].velocity * state[i].mass / SOLAR_MASS;
    }
}

double energy(const body state[BODIES_COUNT]) {
    double energy = 0;

    for (unsigned int i = 0; i < BODIES_COUNT; ++i) {
        const body &body1 = state[i];
        energy += 0.5 * body1.mass * body1.velocity.norm();
        for (unsigned int j = i + 1; j < BODIES_COUNT; ++j) {
            const body &body2 = state[j];
            vector3d r12 = body1.position - body2.position;
            energy -= body1.mass * body2.mass / sqrt(r12.norm());
        }
    }

    return energy;
}

body state[] = {
        // Sun
        {
                .name = "sun",
                .position = {
                        0,
                        0,
                        0
                },
                .velocity = {
                        0,
                        0,
                        0
                },
                .mass = SOLAR_MASS
        },
        // Jupiter
        {
                .name = "jupiter",
                .position = {
                        4.84143144246472090e+00,
                        -1.16032004402742839e+00,
                        -1.03622044471123109e-01
                },
                .velocity = {
                        1.66007664274403694e-03 * DAYS_PER_YEAR,
                        7.69901118419740425e-03 * DAYS_PER_YEAR,
                        -6.90460016972063023e-05 * DAYS_PER_YEAR
                },
                .mass = 9.54791938424326609e-04 * SOLAR_MASS
        },
        // Saturn
        {
                .name = "saturn",
                .position = {
                        8.34336671824457987e+00,
                        4.12479856412430479e+00,
                        -4.03523417114321381e-01
                },
                .velocity = {
                        -2.76742510726862411e-03 * DAYS_PER_YEAR,
                        4.99852801234917238e-03 * DAYS_PER_YEAR,
                        2.30417297573763929e-05 * DAYS_PER_YEAR
                },
                .mass = 2.85885980666130812e-04 * SOLAR_MASS
        },
        // Uranus
        {
                .name = "uranus",
                .position = {
                        1.28943695621391310e+01,
                        -1.51111514016986312e+01,
                        -2.23307578892655734e-01
                },
                .velocity = {
                        2.96460137564761618e-03 * DAYS_PER_YEAR,
                        2.37847173959480950e-03 * DAYS_PER_YEAR,
                        -2.96589568540237556e-05 * DAYS_PER_YEAR
                },
                .mass = 4.36624404335156298e-05 * SOLAR_MASS
        },
        // Neptune
        {
                .name = "neptune",
                .position = {
                        1.53796971148509165e+01,
                        -2.59193146099879641e+01,
                        1.79258772950371181e-01
                },
                .velocity = {
                        2.68067772490389322e-03 * DAYS_PER_YEAR,
                        1.62824170038242295e-03 * DAYS_PER_YEAR,
                        -9.51592254519715870e-05 * DAYS_PER_YEAR
                },
                .mass = 5.15138902046611451e-05 * SOLAR_MASS
        }
};


int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "This is " << argv[0] << std::endl;
        std::cout << "Call this program with an integer as program argument" << std::endl;
        std::cout << "(to set the number of iterations for the n-body simulation)." << std::endl;
        return EXIT_FAILURE;
    } else {
        const unsigned int n = atoi(argv[1]);
        clock_t start_time, end_time;
        start_time = clock();

        offset_momentum(state);
        std::cout << energy(state) << std::endl;

 // uncomment following lines to create csv file
 /*
        std::ofstream outFile;
        outFile.open("body_position_cpp_500000.csv", std::ios::out);
        outFile << "name of the body" << "," << "position x" << "," <<
        "position y" << "," << "position z" << std::endl;
*/
        for (int i = 0; i < n; ++i) {
            advance(state, 0.01);
        }

 // uncomment following lines to create csv file
/*
            for (int j = 0; j < BODIES_COUNT; ++j) {
                outFile << state[j].name << "," << state[j].position.x << "," <<
                state[j].position.y << "," << state[j].position.z << std::endl;
            }
*/
        std::cout << energy(state) << std::endl;
        end_time = clock();
        double duration;
        duration = double(end_time - start_time) / (CLOCKS_PER_SEC);
        std::cout << duration << std::endl;
        return EXIT_SUCCESS;

    }
}