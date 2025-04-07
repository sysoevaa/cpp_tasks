#include <cmath>
#include <type_traits>
#include <vector>
#include <tuple>
#include <iostream>
#include <cstdarg>

namespace geometry {
  const double epsilon = 1e-7;
  const double pi = acos(-1);

  bool isEqual(double first, double second) {
    return std::abs(first - second) < epsilon;
  }

  bool isZero(double value) {
    return std::abs(value) < epsilon;
  }
}

struct Point {
  double x, y;

  Point(double x, double y) : x(x), y(y) {}

  Point& operator-=(const Point& other) {
    x -= other.x;
    y -= other.y;
    return *this;
  }

  Point& operator+=(const Point& other) {
    x += other.x;
    y += other.y;
    return *this;
  }

  double len() const {
    return std::hypot(x, y);
  }

  void rotate(double angle, const Point& center);
};

bool operator==(const Point& pt1, const Point& pt2) {
  return geometry::isEqual(pt1.x, pt2.x) && geometry::isEqual(pt1.y, pt2.y);
}

bool operator!=(const Point& pt1, const Point& pt2) {
  return !(pt1 == pt2);
}

Point operator+(const Point& pt1, const Point& pt2) {
  auto copy = pt1;
  copy += pt2;
  return copy;
}

Point operator-(const Point& pt1, const Point& pt2) {
  auto copy = pt1;
  copy -= pt2;
  return copy;
}

Point operator*(const Point& pt, double k) {
  return {pt.x * k, pt.y * k};
}

Point operator/(const Point& pt, double k) {
  return {pt.x / k, pt.y / k};
}

void Point::rotate(double angle, const Point& center) {
  double prev_x = x;
  double prev_y = y;
  x =
      (prev_x - center.x) * cos(angle) - (prev_y - center.y) * sin(angle) +
      center.x;
  y =
      (prev_x - center.x) * sin(angle) + (prev_y - center.y) * cos(angle) +
      center.y;
}

namespace geometry {
  double dotProduct(const Point& first, const Point& second) {
    return first.x * second.x + first.y * second.y;
  }

  double crossProduct(const Point& first, const Point& second) {
    return first.x * second.y - first.y * second.x;
  }

  double angle(const Point& first, const Point& second) {
    return atan2(crossProduct(first, second), dotProduct(first, second));
  }
}

class Line {
 public:
  Line(double coefficient, double shift) : A(coefficient), B(-1),
      C(shift) {}

  Line(const Point& pt1, const Point& pt2);

  Line(const Point& pt, double coefficient) : A(coefficient), B(-1), C(
      pt.y - coefficient * pt.x) {}

  std::tuple<double, double, double> getCoefficient() const {
    return {A, B, C};
  }

  Point intersection(const Line& other);

 private:
  double A, B, C;
};

Line::Line(const Point& pt1, const Point& pt2) {
  if (pt1.x == pt2.x) {
    A = 1;
    B = 0;
    C = -pt1.x;
    return;
  }
  if (pt1.y == pt2.y) {
    A = 0;
    B = 1;
    C = -pt1.y;
    return;
  }
  A = (pt2.y - pt1.y) / (pt2.x - pt1.x);
  B = -1;
  C = pt2.y - A * pt2.x;
}

Point Line::intersection(const Line& other) {
  auto [A2, B2, C2] = other.getCoefficient();
  double det = A * B2 - B * A2;
  double det1 = -C * B2 + C2 * B;
  double det2 = -A * C2 + A2 * C;
  return Point(det1 / det, det2 / det);
}

bool operator==(const Line& first, const Line& second) {
  auto [A1, B1, C1] = first.getCoefficient();
  auto [A2, B2, C2] = second.getCoefficient();
  return geometry::isEqual(A1 * B2, B1 * A2) && geometry::isEqual(C1, C2);
}

class Shape {
 public:
  virtual void rotate(const Point& center, double angle) = 0;

  virtual void reflect(const Point& center) = 0;

  virtual void reflect(const Line& axis) = 0;

  virtual void scale(const Point& center, double coefficient) = 0;

  virtual double perimeter() const = 0;

  virtual double area() const = 0;

  virtual bool isCongruentTo(const Shape& another) const = 0;

  virtual bool isSimilarTo(const Shape& another) const = 0;

  virtual bool containsPoint(const Point& point) const = 0;

  virtual ~Shape() = default;
};

bool operator==(const Shape& first, const Shape& second) {
  return first.isSimilarTo(second);
}

bool operator!=(const Shape& first, const Shape& second) {
  return !(first == second);
}

class Polygon : public Shape {
 public:
  Polygon(const std::vector<Point>& points) : vertices_(points) {
  }

  Polygon(const auto& ...tail) : vertices_({tail...}) {
  }

  void rotate(const Point& center, double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, double coefficient) override;

  double perimeter() const override;

  double area() const override;

  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  size_t verticesCount() const { return vertices_.size(); }

  const std::vector<Point> getVertices() const { return vertices_; }

  virtual bool isConvex() const;

  ~Polygon() override = default;

 protected:
  Polygon() {}

  static void sort(std::vector<Point>& points);

  std::pair<bool, double> isSimilar(const Shape& another) const;

  std::vector<Point> vertices_;
};

bool operator==(const Polygon& first, const Polygon& second) {
  if (first.verticesCount() != second.verticesCount()) {
    return false;
  }
  auto vert_1 = first.getVertices(), vert_2 = second.getVertices();
  auto cmp = [&](const Point& a, const Point& b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  };
  std::sort(vert_1.begin(), vert_1.end(), cmp);
  std::sort(vert_2.begin(), vert_2.end(), cmp);
  for (size_t i = 0; i < vert_1.size(); ++i) {
    if (vert_1[i] != vert_2[i]) {
      return false;
    }
  }
  return true;
}

bool operator!=(const Polygon& first, const Polygon& second) {
  return !(first == second);
}

void Polygon::rotate(const Point& center, double angle) {
  for (auto& point: vertices_) {
    point.rotate(angle, center);
  }
}

void Polygon::reflect(const Point& center) {
  scale(center, 1);
}

void Polygon::reflect(const Line& axis) {
  for (auto& point: vertices_) {
    auto [A, B, C] = axis.getCoefficient();
    if (B == 0) {
      point.x = point.x + 2 * (-C - point.x);
      continue;
    }
    if (A != 0) {
      point.rotate(atan2(-A, B), {0, 0});
    }
    point.y = point.y + 2 * (-C - point.y);
    if (A != 0) {
      point.rotate(atan2(A, B), {0, 0});
    }
  }
}

void Polygon::scale(const Point& center, double coefficient) {
  for (auto& point: vertices_) {
    Point p(center.x - point.x, center.y - point.y);
    p = p * (coefficient - 1);
    point -= p;
  }
}

double Polygon::perimeter() const {
  double answer = 0;
  for (size_t i = 0; i + 1 < vertices_.size(); ++i) {
    answer += (vertices_[i + 1] - vertices_[i]).len();
  }
  answer += (vertices_.front() - vertices_.back()).len();
  return answer;
}

double Polygon::area() const {
  double answer = 0;
  for (size_t i = 0; i + 1 < vertices_.size(); ++i) {
    answer += geometry::crossProduct(vertices_[i + 1], vertices_[i]);
  }
  answer += geometry::crossProduct(vertices_.front(), vertices_.back());
  answer /= 2;
  return std::abs(answer);
}

bool Polygon::containsPoint(const Point& point) const {
  double angle = 0;
  for (size_t i = 0; i + 1 < vertices_.size(); ++i) {
    auto first = vertices_[i] - point, second = vertices_[i + 1] - point;
    angle += geometry::angle(first, second);
  }
  auto first = vertices_.back() - point, second = vertices_.front() - point;
  angle += geometry::angle(first, second);
  return !geometry::isZero(angle);
}

bool Polygon::isConvex() const {
  if (vertices_.size() <= 3) {
    return true;
  }
  auto copy = vertices_;
  sort(copy);
  auto mod = [&](size_t number) {
    return (number >= verticesCount() ? number - verticesCount() : number);
  };
  for (size_t i = 0; i < copy.size(); ++i) {
    auto first = copy[i] - copy[mod(i + 1)], second =
        copy[mod(i + 2)] - copy[mod(i + 1)];
    if (geometry::crossProduct(second, first) < geometry::epsilon) {
      return false;
    }
  }
  return true;
}

bool Polygon::isSimilarTo(const Shape& another) const {
  auto [flag, attitude] = isSimilar(another);
  return flag;
}

bool Polygon::isCongruentTo(const Shape& another) const {
  auto [flag, attitude] = isSimilar(another);
  return std::abs(attitude - 1) < geometry::epsilon;
}

void Polygon::sort(std::vector<Point>& points) {
  std::sort(points.begin(), points.end(),
            [](const Point& a, const Point& b) {
              return std::make_pair(a.x, a.y) < std::make_pair(b.x, b.y);
            });
  Point point = points[0];
  std::stable_sort(points.begin() + 1, points.end(),
                   [&](const Point& first, const Point& second) {
                     return geometry::crossProduct(first - point,
                                                   (second - point)) > 0;
                   });
}

std::pair<bool, double> Polygon::isSimilar(const Shape& another) const {
  const Polygon* casted = dynamic_cast<const Polygon*>(&another);
  if (casted == nullptr) {
    return {false, -1};
  }
  if (vertices_.size() != casted->vertices_.size()) {
    return {false, -1};
  }
  auto mod_add = [&](size_t number, size_t value) {
    return (number + value >= vertices_.size() ?
            number + value - vertices_.size() : number + value);
  };

  auto vertex = vertices_;
  for (size_t start = 0; start < vertices_.size(); ++start) {
    bool ok = true;
    double attitude = 0;
    for (size_t i = 0; i < vertices_.size(); ++i) {
      Point first = vertices_[mod_add(i, start + 1)] -
                    vertices_[mod_add(i, start)];
      Point second =
          casted->vertices_[mod_add(i, 1)] - casted->vertices_[i];
      double first_len = first.len();
      double second_len = second.len();
      if (attitude == 0) {
        attitude = second_len / first_len;
      } else {
        if (!geometry::isEqual(attitude, second_len / first_len)) {
          ok = false;
          break;
        }
      }
    }
    if (ok) {
      return {true, attitude};
    }
  }
  std::reverse(vertex.begin(), vertex.end());
  for (size_t start = 0; start < vertex.size(); ++start) {
    bool ok = true;
    double attitude = 0;
    for (size_t i = 0; i < vertex.size(); ++i) {
      Point first = vertex[mod_add(i, start + 1)] -
                    vertex[mod_add(i, start)];
      Point second =
          casted->vertices_[mod_add(i, 1)] - casted->vertices_[i];
      double first_len = first.len();
      double second_len = second.len();
      if (attitude == 0) {
        attitude = second_len / first_len;
      } else {
        if (std::abs(attitude - second_len / first_len) > 0.00001) {
          ok = false;
          break;
        }
      }
    }
    if (ok) {
      return {true, attitude};
    }
  }
  return {false, -1};
}

class Ellipse : public Shape {
 public:
  Ellipse(const Point& focus_1, const Point& focus_2, double dist) : focuses_{
      focus_1, focus_2}, eccentricity_((focus_2 - focus_1).len() / dist),
      semi_major_axis_(dist / 2) {
  }

  void rotate(const Point& center, double angle) override;

  void reflect(const Point& center) override;

  void reflect(const Line& axis) override;

  void scale(const Point& center, double coefficient) override;

  double perimeter() const override;

  double area() const override;

  double semiMajorAxis() const { return semi_major_axis_; }

  bool isCongruentTo(const Shape& another) const override;

  bool isSimilarTo(const Shape& another) const override;

  bool containsPoint(const Point& point) const override;

  std::pair<Line, Line> directrices() const;

  std::pair<Point, Point> focuses() const;

  double eccentricity() const { return eccentricity_; }

  Point center() const { return (focuses_[0] + focuses_[1]) / 2; }

  ~Ellipse() override = default;

 protected:
  Point focuses_[2];
  double eccentricity_;
  double semi_major_axis_;
};

void Ellipse::rotate(const Point& center, double angle) {
  for (auto& focuse: focuses_) {
    focuse.rotate(angle, center);
  }
}

void Ellipse::reflect(const Point& center) {
  scale(center, 1);
}

void Ellipse::reflect(const Line& axis) {
  Polygon p(focuses_[0], focuses_[1]);
  p.reflect(axis);
  auto rotated = p.getVertices();
  focuses_[0] = rotated[0];
  focuses_[1] = rotated[1];
}

void Ellipse::scale(const Point& center, double coefficient) {
  Polygon p(focuses_[0], focuses_[1]);
  p.scale(center, coefficient);
  auto rotated = p.getVertices();
  focuses_[0] = rotated[0];
  focuses_[1] = rotated[1];
  semi_major_axis_ *= coefficient;
}

double Ellipse::perimeter() const {
  return std::comp_ellint_2(eccentricity_) * 4 * semi_major_axis_;
  // return 0;
}

double Ellipse::area() const {
  return semi_major_axis_ * semi_major_axis_ * geometry::pi *
         sqrt(1 - eccentricity_ * eccentricity_);
}

bool Ellipse::isCongruentTo(const Shape& another) const {
  auto casted = dynamic_cast<const Ellipse*>(&another);
  if (casted == nullptr) {
    return false;
  }
  return semi_major_axis_ == casted->semi_major_axis_ &&
         eccentricity_ == casted->eccentricity_;
}

bool Ellipse::isSimilarTo(const Shape& another) const {
  auto casted = dynamic_cast<const Ellipse*>(&another);
  if (casted == nullptr) {
    return false;
  }
  return eccentricity_ == casted->eccentricity_;
}

bool Ellipse::containsPoint(const Point& point) const {
  auto val = (point - focuses_[0]).len() + (point - focuses_[1]).len();
  return val < 2 * semi_major_axis_;
}

std::pair<Line, Line> Ellipse::directrices() const {
  return {Line(Point(semi_major_axis_ / eccentricity_, 1),
               Point(semi_major_axis_ / eccentricity_, 2)),
          Line(Point(-semi_major_axis_ / eccentricity_, 1),
               Point(-semi_major_axis_ / eccentricity_, 2))};
}

std::pair<Point, Point> Ellipse::focuses() const {
  return {focuses_[0], focuses_[1]};
}

bool operator==(const Ellipse& first, const Ellipse& second) {
  return first.focuses() == second.focuses() &&
         first.semiMajorAxis() == second.semiMajorAxis();
}

bool operator!=(const Ellipse& first, const Ellipse& second) {
  return !(first == second);
}

class Circle : public Ellipse {
 public:
  Circle(const Point& center, const double radius) : Ellipse(center, center,
                                                             2 * radius) {
  }

  double radius() const { return semi_major_axis_; }

  ~Circle() override = default;

 private:
};

class Rectangle : public Polygon {
 public:
  Rectangle(const Point& vert_1, const Point& vert_2, double coefficient)
      : Polygon() {
    Point first = vert_1, second = vert_2;
    if (first.x > second.x) {
      std::swap(first, second);
    }
    if (coefficient < 1 - geometry::epsilon) {
      coefficient = 1 / coefficient;
    }
    Point diagonal = first - second;
    double less = diagonal.len() / std::hypot(1, coefficient);
    diagonal.x /= diagonal.len();
    diagonal.y /= diagonal.len();
    diagonal.rotate(atan2(coefficient, 1), second);
    diagonal.x *= less * coefficient;
    diagonal.y *= less * coefficient;
    vertices_.push_back(second + diagonal);
    vertices_.push_back(first - diagonal);
    vertices_.push_back(vert_1);
    vertices_.push_back(vert_2);
    Polygon::sort(vertices_);
  }

  Point center() const;

  std::pair<Line, Line> diagonals() const;

  bool isConvex() const override { return true; }

  double area() const override;

  ~Rectangle() override = default;
};

double Rectangle::area() const {
  auto edge1 = vertices_[1] - vertices_[0];
  auto edge2 = vertices_[2] - vertices_[1];
  return edge1.len() * edge2.len();
}

Point Rectangle::center() const {
  return Point((vertices_[0].x + vertices_[2].x) / 2,
               (vertices_[0].y + vertices_[2].y) / 2);
}

std::pair<Line, Line> Rectangle::diagonals() const {
  return {Line(vertices_[0], vertices_[2]), Line(vertices_[1], vertices_[3])};
}

class Square : public Rectangle {
 public:
  Square(const Point& vert_1, const Point& vert_2) : Rectangle(vert_1, vert_2,
                                                               1) {
  }

  Circle circumscribedCircle() const;

  Circle inscribedCircle() const;

  ~Square() override = default;
};

Circle Square::circumscribedCircle() const {
  Point first = vertices_[0], second = vertices_[2];
  return {Point(first.x + (second.x - first.x) / 2,
                first.y + (second.y - first.y) / 2),
          (second - first).len()};
}

Circle Square::inscribedCircle() const {
  Point first = vertices_[0], second = vertices_[2];
  return {Point(first.x + (second.x - first.x) / 2,
                first.y + (second.y - first.y) / 2),
          second.x - first.x};
}

class Triangle : public Polygon {
 public:
  Triangle(const Point& vert_1, const Point& vert_2, const Point& vert_3)
      : Polygon(vert_1, vert_2, vert_3) {
  }

  Circle circumscribedCircle() const;

  Circle inscribedCircle() const;

  Point centroid() const;

  Point orthocenter() const;

  Line EulerLine() const;

  Circle ninePointsCircle() const;

  bool isConvex() const override { return true; }

  double area() const override;

  ~Triangle() override = default;
};

double Triangle::area() const {
  auto edge1 = vertices_[1] - vertices_[0], edge2 = vertices_[2] - vertices_[1];
  return std::abs(geometry::crossProduct(edge1, edge2)) / 2;
}

Circle Triangle::inscribedCircle() const {
  Point A = vertices_[0], B = vertices_[1], C = vertices_[2];
  double a = (B - C).len(), b = (A - C).len(), c = (A - B).len();
  Point Ap = C + (B - C) / (B - C).len() * a * b / (c + b);
  Point Cp = A + (B - A) / (B - A).len() * b * c / (a + b);
  Line A_biss(A, Ap), C_biss(C, Cp);
  auto center = A_biss.intersection(C_biss);
  return Circle(center, area() / Polygon::perimeter() * 2);

}

Circle Triangle::circumscribedCircle() const {
  Point A = vertices_[0], B = vertices_[1], C = vertices_[2];
  double R = (A - C).len() / 2 /
             sin(atan2(std::abs(geometry::crossProduct(A - B, B - C)),
                       std::abs(geometry::dotProduct(A - B, B - C))));
  Point AB_center((A.x + B.x) / 2, (A.y + B.y) / 2);
  double center_to_half = sqrt(
      R * R - (A - AB_center).len() * (A - AB_center).len());
  Line AB(A, B);
  auto [a, b, c] = AB.getCoefficient();
  Point perp(a / sqrt(a * a + b * b), b / sqrt(a * a + b * b));
  Circle c1(AB_center + perp * center_to_half, R);
  Circle c2(AB_center - perp * center_to_half, R);
  auto [xc1, yc1] = c1.center();
  if (((C.x - xc1) * (C.x - xc1) + (C.y - yc1) * (C.y - yc1) -
       c1.radius() * c1.radius()) < geometry::epsilon) {
    return c1;
  }
  return Circle(AB_center - perp * center_to_half, R);
}

Point Triangle::centroid() const {
  Line first(vertices_[0], Point((vertices_[1].x + vertices_[2].x) / 2,
                                 (vertices_[1].y + vertices_[2].y) / 2));
  Line second(vertices_[2], Point((vertices_[0].x + vertices_[1].x) / 2,
                                  (vertices_[0].y + vertices_[1].y) / 2));
  return first.intersection(second);
}

Point Triangle::orthocenter() const {
  Point A = vertices_[0], B = vertices_[1], C = vertices_[2];
  double angle_a = geometry::angle(C - A, B - A);
  double angle_c = geometry::angle(A - C, B - C);
  Point Ap = C + (B - C) / (B - C).len() * (A - C).len() * cos(angle_c);
  Point Cp = A + (B - A) / (B - A).len() * (A - C).len() * cos(angle_a);
  Line first(C, Cp), second(A, Ap);
  return first.intersection(second);
}

Line Triangle::EulerLine() const {
  return Line(centroid(), ninePointsCircle().center());
}

Circle Triangle::ninePointsCircle() const {
  Point A = vertices_[0], B = vertices_[1], C = vertices_[2];
  Point center1 = {(A.x + B.x) / 2, (A.y + B.y) / 2},
      center2 = {(A.x + C.x) / 2, (A.y + C.y) / 2}, center3 = {(B.x + C.x) / 2,
                                                               (B.y + C.y) / 2};
  return Triangle(center1, center2, center3).circumscribedCircle();
}
