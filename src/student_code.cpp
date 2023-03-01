#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    // TODO Part 1.
    std::vector<Vector2D> next_level;
    for (size_t i = 0; i < points.size() - 1; ++i) {
        Vector2D new_point((1 - t) * points[i] + t * points[i+1]);
        next_level.push_back(new_point);
    }
    return next_level;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> next_level;
      for (size_t i = 0; i < points.size() - 1; ++i) {
          Vector3D new_point((1 - t) * points[i] + t * points[i + 1]);
          next_level.push_back(new_point);
      }
      return next_level;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    // TODO Part 2.
      std::vector<Vector3D> level;
    for (int i = 0; i < points.size() - 1; i++) {
        if (i == 0) {
            level = evaluateStep(points, t);
        }
        else {
            level = evaluateStep(level, t);
        }
    }
    return level[0];
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    // TODO Part 2.
      std::vector<Vector3D> moving_curve;
      for (size_t i = 0; i < controlPoints.size(); ++i) {
          Vector3D new_point(evaluate1D(controlPoints[i], u));
          moving_curve.push_back(new_point);
      }
      return evaluate1D(moving_curve, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.
      HalfedgeCIter e(this->halfedge());
      Vector3D A(position);
      Vector3D norm(0, 0, 0);
      while (true) {
          e = e->next();
          Vector3D B(e->vertex()->position);
          e = e->next();
          Vector3D C(e->vertex()->position);
          norm += (cross(B - A, C - A));
          e = e->twin();
          if (e == (this->halfedge())) { 
              break;
          }
    }
    norm.normalize();
    return norm;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.
    if (e0->isBoundary()) goto finish;
      
    HalfedgeIter BC(e0->halfedge());
    HalfedgeIter CB(BC->twin());

    HalfedgeIter CD(BC->next());
    HalfedgeIter DC(CD->twin());
    HalfedgeIter DB(CD->next());
    HalfedgeIter BD(DB->twin());

    HalfedgeIter BA(CB->next());
    HalfedgeIter AB(BA->twin());
    HalfedgeIter AC(BA->next());
    HalfedgeIter CA(AC->twin());

    EdgeIter A_B(AB->edge());
    EdgeIter B_D(BD->edge());
    EdgeIter C_D(DC->edge());
    EdgeIter A_C(CA->edge());
    EdgeIter B_C(BC->edge());

    FaceIter ABC(CB->face());
    FaceIter BCD(BC->face());

    VertexIter D = DC->vertex();
    VertexIter C = CB->vertex();
    VertexIter B = BC->vertex();
    VertexIter A = AB->vertex();

    BC->setNeighbors(CD, CB, A, B_C, BCD);
    CD->setNeighbors(DB, BD, D, B_D, BCD);
    DB->setNeighbors(BC, AB, B, A_B, BCD);

    CB->setNeighbors(BA, BC, D, B_C, ABC);
    BA->setNeighbors(AC, CA, A, A_C, ABC);
    AC->setNeighbors(CB, DC, C, C_D, ABC);

    CA->setNeighbors(CA->next(), BA, C, A_C, CA->face());
    DC->setNeighbors(DC->next(), AC, D, C_D, DC->face());
    BD->setNeighbors(BD->next(), CD, B, B_D, BD->face());
    AB->setNeighbors(AB->next(), DB, A, A_B, AB->face());

    D->halfedge() = CB;
    C->halfedge() = AC;
    B->halfedge() = DB;
    A->halfedge() = BC;

    A_B->halfedge() = AB;
    A_C->halfedge() = CA;
    B_C->halfedge() = BC;
    B_D->halfedge() = CD;
    C_D->halfedge() = DC;

    ABC->halfedge() = BA;
    BCD->halfedge() = CD;

    finish:
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.
      VertexIter M = newVertex();
      if (e0->isBoundary()) goto finish;
      HalfedgeIter BC(e0->halfedge());
      HalfedgeIter CB(BC->twin());

      HalfedgeIter CD(BC->next());
      HalfedgeIter DC(CD->twin());
      HalfedgeIter DB(CD->next());
      HalfedgeIter BD(DB->twin());

      HalfedgeIter BA(CB->next());
      HalfedgeIter AB(BA->twin());
      HalfedgeIter AC(BA->next());
      HalfedgeIter CA(AC->twin());

      EdgeIter A_B(AB->edge());
      EdgeIter B_D(BD->edge());
      EdgeIter C_D(DC->edge());
      EdgeIter A_C(CA->edge());

      FaceIter ABC(CB->face());
      FaceIter BCD(BC->face());

      VertexIter D = DC->vertex();
      VertexIter C = CB->vertex();
      VertexIter B = BC->vertex();
      VertexIter A = AB->vertex();

      EdgeIter M_A = newEdge();
      EdgeIter M_B = newEdge();
      EdgeIter M_C(BC->edge());
      EdgeIter M_D = newEdge();

      FaceIter MAB = newFace();
      FaceIter MCD = newFace();

      HalfedgeIter MC = newHalfedge();
      HalfedgeIter MB = newHalfedge();
      HalfedgeIter MA = newHalfedge();
      HalfedgeIter MD = newHalfedge();
      HalfedgeIter DM = newHalfedge();
      HalfedgeIter AM = newHalfedge();

      AC->setNeighbors(CB, CA, A, A_C, ABC);
      CB->setNeighbors(MA, MC, C, M_C, ABC);
      MA->setNeighbors(AC, AM, M, M_A, ABC);

      CD->setNeighbors(DM, DC, C, C_D, MCD);
      DM->setNeighbors(MC, MD, D, M_D, MCD);
      MC->setNeighbors(CD, CB, M, M_C, MCD);

      DB->setNeighbors(BC, BD, D, B_D, BCD);
      BC->setNeighbors(MD, MB, B, M_B, BCD);
      MD->setNeighbors(DB, DM, M, M_D, BCD);

      BA->setNeighbors(AM, AB, B, A_B, MAB);
      AM->setNeighbors(MB, MA, A, M_A, MAB);
      MB->setNeighbors(BA, BC, M, M_B, MAB);

      D->halfedge() = DB;
      C->halfedge() = CD;
      B->halfedge() = BA;
      A->halfedge() = AC;
      M->halfedge() = MA;

      A_B->halfedge() = AB;
      A_C->halfedge() = CA;
      B_D->halfedge() = BD;
      C_D->halfedge() = DC;

      M_A->halfedge() = MA;
      M_B->halfedge() = MB;
      M_C->halfedge() = MC;
      M_D->halfedge() = MD;

      M_A->isNew = true;
      M_D->isNew = true;

      M->isNew = true;

      ABC->halfedge() = CB;
      BCD->halfedge() = BC;
      MAB->halfedge() = AM;
      MCD->halfedge() = DM;

      M->position = (B->position + C->position) / 2;

    finish:
    return M;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    
    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    
    // 4. Flip any new edge that connects an old and new vertex.

    // 5. Copy the new vertex positions into final Vertex::position.
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          HalfedgeIter e = v->halfedge();
          Vector3D neighbor_sum(0, 0, 0);
          int n = 0;
          v->isNew = false;
          while(true) {
              neighbor_sum += e->next()->vertex()->position;
              e = e->next()->next()->twin();
              n++;
              if (e == v->halfedge()) {
                  break;
              }
          }
          float u = ((n == 3) ? (static_cast<float>(3) / 16) : (3 / (static_cast<float>(8) * n)));
          v->newPosition = (1 - n * u) * v->position + u * neighbor_sum;
      }
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          Vector3D new_position;
          HalfedgeIter AB = e->halfedge();
          HalfedgeIter BA = AB->twin();
          HalfedgeIter CA = AB->next()->next();
          HalfedgeIter DB = BA->next()->next();

          VertexIter A = AB->vertex();
          VertexIter B = BA->vertex();
          VertexIter C = CA->vertex();
          VertexIter D = DB->vertex();

          new_position = (static_cast<double>(3) / 8) * (A->position + B->position) + (static_cast<double>(1) / 8) * (C->position + D->position);

          e->newPosition = new_position;
      }
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if (!((e->halfedge()->vertex()->isNew) || (e->halfedge()->twin()->vertex()->isNew))) {
              Vector3D position(e->newPosition);
              VertexIter M = mesh.splitEdge(e);
              M->newPosition = position;
          }
      }
      for (EdgeIter e = mesh.edgesBegin(); e != mesh.edgesEnd(); e++) {
          if (e->isNew) {
              if ((e->halfedge()->vertex()->isNew) != (e->halfedge()->twin()->vertex()->isNew)) {
                  mesh.flipEdge(e);
              }
          }
          e->isNew = false;
      }
      for (VertexIter v = mesh.verticesBegin(); v != mesh.verticesEnd(); v++) {
          v->position = v->newPosition;
          v->isNew = false;
      }
  }
}
