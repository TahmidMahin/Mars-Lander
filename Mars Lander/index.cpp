#include<bits/stdc++.h>

using namespace std;

const float PI = 3.14159265358979323846;

class Point {
public:
    float x, y;
    Point() : x(-1), y(-1) {}
    Point(float x, float y) : x(x), y(y) {}
    float distance(Point p) {
        return sqrt((x - p.x) * (x - p.x) + (y - p.y) * (y - p.y));
    }
};

class LineSegment {
public:
    Point p1, p2;
    LineSegment() : p1(), p2() {}
    LineSegment(Point p1, Point p2) : p1(p1), p2(p2) {}
    float length() {
        return p1.distance(p2);
    }
    Point checkCollision(LineSegment l2) {
        float x1 = p1.x, y1 = p1.y, x2 = p2.x, y2 = p2.y;
        float x3 = l2.p1.x, y3 = l2.p1.y, x4 = l2.p2.x, y4 = l2.p2.y;
        float den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4);
        if (den == 0) {
            return Point(-1, -1);
        }
        float t = ((x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)) / den;
        float u = -((x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)) / den;
        if (t >= 0 && t <= 1 && u >= 0 && u <= 1) {
            return Point(x1 + t * (x2 - x1), y1 + t * (y2 - y1));
        }
        return Point(-1, -1);
    }
};

class Shuttle {
public:
    Point p;
    float vx;
    float vy;
    float fuel;
    float angle;
    float power;
    Shuttle() : p(), vx(0), vy(0), fuel(0), angle(0), power(0) {}
    Shuttle(Point p, float vx, float vy, float fuel, float angle, float power) : p(p), vx(vx), vy(vy), fuel(fuel), angle(angle), power(power) {}
    void move(float angle, int thrust) {
        if (angle > this->angle) {
            if (angle - this->angle > 15)
                this->angle += 15;
            else
                this->angle = angle;
        } else if (angle < this->angle) {
            if(this->angle - angle > 15)
                this->angle -= 15;
            else
                this->angle = angle;
        }

        if (thrust > this->power)
            power += min(thrust - this->power, 1.0f);
        else if (thrust < this->power)
            power -= min(this->power - thrust, 1.0f);

        fuel -= power;

        float ax = power * sin(angle * PI / 180);
        float ay = power * cos(angle * PI / 180) - 3.771;

        vx += ax;
        vy += ay;

        p.x += vx;
        p.y += vy;
    }
};

class Matrix {
public:
    int row, column;
    float **data;
    Matrix() : row(0), column(0) {}
    Matrix(int row, int column) : row(row), column(column) {
        data = new float*[row];
        for (int i = 0; i < row; i++) {
            data[i] = new float[column];
        }
    }
    void set(int i, int j, float value) {
        data[i][j] = value;
    }
    float get(int i, int j) {
        return data[i][j];
    }
    void print() {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                cout << data[i][j] << " ";
            }
            cout << endl;
        }
    }
    Matrix operator*(Matrix m) {
        Matrix result(row, m.column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < m.column; j++) {
                float sum = 0;
                for (int k = 0; k < column; k++) {
                    sum += data[i][k] * m.data[k][j];
                }
                result.set(i, j, sum);
            }
        }
        return result;
    }
    Matrix operator*(float f) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] * f);
            }
        }
        return result;
    }
    Matrix operator/(float f) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] / f);
            }
        }
        return result;
    }
    Matrix operator+(Matrix m) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] + m.data[i][j]);
            }
        }
        return result;
    }
    Matrix operator-(Matrix m) {
        Matrix result(row, column);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(i, j, data[i][j] - m.data[i][j]);
            }
        }
        return result;
    }
    Matrix transpose() {
        Matrix result(column, row);
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < column; j++) {
                result.set(j, i, data[i][j]);
            }
        }
        return result;
    }

};

vector<LineSegment> surface;
LineSegment flatGround;

float sigmoid(float x) {
    return 1 / (1 + exp(-x));
}

class Model {
public:
    Matrix W, b;
    vector<Matrix> cacheX;
    vector<Matrix> cacheA;
    vector<Matrix> gradX;

    Model() {
        W = Matrix(2, 7);
        b = Matrix(2, 1);
    }

    void initialize() {
        for (int i=0; i<2; i++) {
            for (int j=0; j<7; j++) {
                W.set(i, j, 2*(rand()/(float)RAND_MAX)-1.0);
            }
        }
        for (int i=0; i<2; i++) {
            b.set(i, 0, 0);
        }
    }

    Matrix makeVector(Shuttle shuttle) {
        Matrix X(7, 1);
        X.set(0, 0, shuttle.p.x/3500-1);
        X.set(1, 0, shuttle.p.y/1500-1);
        X.set(2, 0, shuttle.vx/500);
        X.set(3, 0, shuttle.vy/500);
        X.set(4, 0, shuttle.fuel/1000-1);
        X.set(5, 0, shuttle.angle/90);
        X.set(6, 0, shuttle.power/4);
        return X;
    } 

    Shuttle makeShuttle(Matrix X) {
        Shuttle s;
        s.p.x = X.get(0, 0) * 3500 + 1;
        s.p.y = X.get(1, 0) * 1500 + 1;
        s.vx = X.get(2, 0) * 500;
        s.vy = X.get(3, 0) * 500;
        s.fuel = X.get(4, 0) * 1000 + 1;
        s.angle = X.get(5, 0) * 90;
        s.power = X.get(6, 0) * 4;
        return s;
    }

    void updateParameter(Matrix deltaW, Matrix deltab, float lr) {
        W = W - deltaW * lr;
        b = b - deltab * lr;
    }

    float lossFunction(Shuttle shuttle, Point p) {
        float loss = p.distance(flatGround.p1) + p.distance(flatGround.p2) - flatGround.length();
        loss += max(abs(shuttle.vx) - 20, 0.0f);
        loss += max(abs(shuttle.vy) - 40, 0.0f);
        loss += tan(abs(shuttle.angle) * PI / 180);
        return loss;
    }


    void forwardPass(Shuttle s) {
        int iteration = 0;

        while (true) {
            iteration++;
            Shuttle s0 = s;
            Matrix X = makeVector(s);

            cacheX.push_back(X);

            Matrix Y = W * X + b;

            float a0 = tanh(Y.get(0, 0));
            float a1 = sigmoid(Y.get(1, 0));
            Matrix A(2, 1);
            A.set(0, 0, a0);
            A.set(1, 0, a1);

            cacheA.push_back(A);
        
            float angle = 90*a0;
            float thrust = floor(5*a1);

            s.move(angle, thrust);

            bool collide = false;
            Point p;
            for (LineSegment l : surface) {
                p = l.checkCollision(LineSegment(s0.p, s.p));
                if (p.x != -1 and p.y != -1) { // add approximation
                    collide = true;
                }
            }

            if (collide) {
                break;
            }
            else {
                X = makeVector(s);
            }
        }
    }

    void backwardPass(float learningRate) {
        int iteration = cacheX.size() - 1;
        Shuttle s0 = makeShuttle(cacheX[iteration]);
        Shuttle s1 = s0, s2 = s0, s3 = s0;
        Point p, q;
        float a0 = cacheA[iteration].get(0, 0);
        float a1 = cacheA[iteration].get(1, 0);
        
        float angle = 90*a0;
        float thrust = floor(5*a1);
        s1.move(angle, thrust);
        for (LineSegment l : surface) {
            p = l.checkCollision(LineSegment(s0.p, s1.p));
            if (p.x != -1 and p.y != -1) { // add approximation
                q = p;
            }
        }
        float loss = lossFunction(s1, q);
        cout << "Loss: " << loss << endl;

        float angle1 = 90*(a0 + 0.1); // derivative of angle
        s2.move(angle1, thrust);
        for (LineSegment l : surface) {
            p = l.checkCollision(LineSegment(s0.p, s2.p));
            if (p.x != -1 and p.y != -1) { // add approximation
                q = p;
            }
        }
        float da0 = (lossFunction(s2, q) - loss) / 0.1;

        float thrust1 = floor(5*(a1 + 0.1)); // derivative of thrust
        s3.move(angle, thrust1);
        for (LineSegment l : surface) {
            p = l.checkCollision(LineSegment(s0.p, s3.p));
            if (p.x != -1 and p.y != -1) { // add approximation
                q = p;
            }
        }
        float da1 = (lossFunction(s3, q) - loss) / 0.1;

        Matrix dY(2, 1);
        dY.set(0, 0, da0*(1-a0*a0));
        dY.set(1, 0, da1*a1*(1-a1));
        Matrix dW = dY * cacheX[iteration].transpose();
        Matrix db = dY;
        updateParameter(dW, db, learningRate);
        Matrix dX = W.transpose() * dY;
        
        for (int i = iteration - 1; i >= 0; i--) {
            s0 = makeShuttle(cacheX[i]);
            s1 = s0;
            s2 = s0;
            s3 = s0;
            a0 = cacheA[i].get(0, 0);
            a1 = cacheA[i].get(1, 0);
            angle = 90*a0;
            thrust = floor(5*a1);
            s1.move(angle, thrust);
            Matrix X = makeVector(s1);

            float angle1 = 90*(a0 + 0.1); // derivative of angle
            s2.move(angle1, thrust);
            Matrix X0 = makeVector(s2);

            float thrust1 = floor(5*(a1 + 0.1)); // derivative of thrust
            s3.move(angle, thrust1);
            Matrix X1 = makeVector(s3);
            
            Matrix dA0 = dX.transpose()*((X0 - X) / 0.1);
            Matrix dA1 = dX.transpose()*((X1 - X) / 0.1);

            dY.set(0, 0, dA0.get(0, 0)*(1-a0*a0));
            dY.set(1, 0, dA1.get(0, 0)*a1*(1-a1));

            dW = dY * cacheX[i].transpose();
            db = dY;
            updateParameter(dW, db, learningRate);
            dX = W.transpose() * dY;

            
        }
    }

    void fit(Shuttle s, int epochs, float learningRate) {
        initialize();
        for (int i = 0; i < epochs; i++) {
            cacheX.clear();
            cacheA.clear();
            forwardPass(s);
            backwardPass(learningRate);
        }
    }
};

int main() {
    freopen("input.txt", "r", stdin);
    int n;
    cin >> n;
    vector<Point> points;
    for (int i = 0; i < n; i++) {
        int x, y;
        cin >> x >> y;
        points.push_back(Point(x, y));
    }

    for (int i = 0; i < n-1; i++) {
        surface.push_back(LineSegment(points[i], points[i+1]));
        if (points[i].y == points[i+1].y) {
            flatGround = LineSegment(points[i], points[i+1]);
        }
    }

    float x, y, vx, vy, fuel, angle, power;

    cin >> x >> y >> vx >> vy >> fuel >> angle >> power;
    Shuttle s = Shuttle(Point(x, y), vx, vy, fuel, angle, power);

    Model model;
    model.fit(s, 200, 0.05);
    
    return 0;
}