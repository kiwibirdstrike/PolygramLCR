#include <Rcpp.h>
using namespace Rcpp;

// functions for computing barycentric coordinates
double triangle_area(NumericMatrix v);
NumericVector rect2bary(NumericVector u, NumericMatrix v);

// functions for linear hat basis
double hat(NumericVector u, NumericMatrix vertex, NumericMatrix star);
NumericMatrix hat_basis_linear(NumericMatrix grid, NumericMatrix vertex, List star);

// [[Rcpp::export]]
double triangle_area(NumericMatrix v)
{
   // v : 3 x 2 matrix
   double d_12, d_13, d_23;
   d_12 = v(0, 0) * v(1, 1) - v(0, 1) * v(1, 0);
   d_13 = v(0, 0) * v(2, 1) - v(0, 1) * v(2, 0);
   d_23 = v(1, 0) * v(2, 1) - v(1, 1) * v(2, 0);
   return 0.5 * (d_23 - d_13 + d_12);
}

// [[Rcpp::export]]
NumericVector rect2bary(NumericVector u, NumericMatrix v)
{
   double area_triangle;
   NumericVector b(3);
   int j;
   //
   area_triangle = triangle_area(v);
   for (j = 0; j < 3; j++)
   {
      NumericMatrix vj = clone(v);
      vj(j, 0) = u[0];
      vj(j, 1) = u[1];
      b[j] = triangle_area(vj) / area_triangle;
   }
   for (int k = 0; k < 3; k++)
      if (std::abs(b[k]) < 1e-10)
         b[k] = 0;
   return b;
}

// [[Rcpp::export]]
double hat(NumericVector u, NumericMatrix vertex, NumericMatrix star)
{
   double s = 0;
   int t;
   NumericMatrix triangle(3, 2);
   NumericVector b(3);
   NumericVector star_row_t(3);
   for (t = 0; t < star.nrow(); t++)
   {
      // triangle
      star_row_t = star(t, _);
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 2; j++)
            triangle(i, j) = vertex(star_row_t[i] - 1, j);
      // IS_IN_TRIANGLE
      b = rect2bary(u, triangle);
      if (b[0] >= 0)
         if (b[1] >= 0)
            if (b[2] >= 0)
            {
               NumericMatrix uv = clone(triangle);
               uv(0, 0) = u[0];
               uv(0, 1) = u[1];
               s = triangle_area(uv) / triangle_area(triangle);
               break;
            }
   }
   return s;
}

// [[Rcpp::export]]
NumericMatrix hat_basis_linear(NumericMatrix grid, NumericMatrix vertex, List star)
{
   int sample_size = grid.nrow();
   int dimension = vertex.nrow();
   NumericMatrix basis(sample_size, dimension);
   int i, j;
   for (i = 0; i < sample_size; i++)
   {
      NumericVector p = grid(i, _);
      for (j = 0; j < dimension; j++)
      {
         NumericMatrix star_j = star[j];
         basis(i, j) = hat(p, vertex, star_j);
      }
   }
   return basis;
}
