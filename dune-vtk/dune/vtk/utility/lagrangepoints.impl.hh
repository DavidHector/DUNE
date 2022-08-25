#pragma once

#include <cassert>
#include <array>

#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/geometry/type.hh>
#include <dune/localfunctions/lagrange/equidistantpoints.hh>

namespace Dune {
namespace Vtk {
namespace Impl {

/**
 *  The implementation of the point set builder is directly derived from VTK.
 *  Modification are a change in data-types and naming scheme. Additionally
 *  a LocalKey is created to reflect the concept of a Dune PointSet.
 *
 *  Included is the license of the BSD 3-clause License included in the VTK
 *  source code from 2020/04/13 in commit b90dad558ce28f6d321420e4a6b17e23f5227a1c
 *  of git repository https://gitlab.kitware.com/vtk/vtk.
 *
    Program:   Visualization Toolkit
    Module:    Copyright.txt

    Copyright (c) 1993-2015 Ken Martin, Will Schroeder, Bill Lorensen
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
      of any contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
    ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/


template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,0>::operator() (GeometryType gt, int /*order*/, Points& points) const
{
  assert(gt.isVertex());
  points.push_back(LP{Vec{},Key{0,0,0}});
}


template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,1>::operator() (GeometryType gt, int order, Points& points) const
{
  assert(gt.isLine());

  // Vertex nodes
  points.push_back(LP{Vec{0.0}, Key{0,dim,0}});
  points.push_back(LP{Vec{1.0}, Key{1,dim,0}});

  if (order > 1) {
    // Inner nodes
    Vec p{0.0};
    for (int i = 0; i < order-1; i++)
    {
      p[0] += 1.0 / order;
      points.push_back(LP{p,Key{0,dim-1,(unsigned int)(i)}});
    }
  }
}


template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,2>::operator() (GeometryType gt, int order, Points& points) const
{
  std::size_t nPoints = numLagrangePoints(gt, order);

  if (gt.isTriangle())
    buildTriangle(nPoints, order, points);
  else if (gt.isQuadrilateral())
    buildQuad(nPoints, order, points);
  else {
    DUNE_THROW(Dune::NotImplemented,
      "Lagrange points not yet implemented for this GeometryType.");
  }

  assert(points.size() == nPoints);
}


// Construct the point set in a triangle element.
// Loop from the outside to the inside
template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,2>::buildTriangle (std::size_t nPoints, int order, Points& points) const
{
  points.reserve(nPoints);

  const int nVertexDOFs = 3;
  const int nEdgeDOFs = 3 * (order-1);

  static const unsigned int vertexPerm[3] = {0,1,2};
  static const unsigned int edgePerm[3]   = {0,2,1};

  auto calcKey = [&](int p) -> Key
  {
    if (p < nVertexDOFs) {
      return Key{vertexPerm[p], dim, 0};
    }
    else if (p < nVertexDOFs+nEdgeDOFs) {
      unsigned int entityIndex = (p - nVertexDOFs) / (order-1);
      unsigned int index = (p - nVertexDOFs) % (order-1);
      return Key{edgePerm[entityIndex], dim-1, index};
    }
    else {
      unsigned int index = p - (nVertexDOFs + nEdgeDOFs);
      return Key{0, dim-2, index};
    }
  };

  std::array<int,3> bindex;

  K order_d = K(order);
  for (std::size_t p = 0; p < nPoints; ++p) {
    barycentricIndex(p, bindex, order);
    Vec point{bindex[0] / order_d, bindex[1] / order_d};
    points.push_back(LP{point, calcKey(p)});
  }
}


// "Barycentric index" is a triplet of integers, each running from 0 to
// <Order>. It is the index of a point on the triangle in barycentric
// coordinates.
template <class K>
void LagrangePointSetBuilder<K,2>::barycentricIndex (int index, std::array<int,3>& bindex, int order)
{
  int max = order;
  int min = 0;

  // scope into the correct triangle
  while (index != 0 && index >= 3 * order)
  {
    index -= 3 * order;
    max -= 2;
    min++;
    order -= 3;
  }

  // vertex DOFs
  if (index < 3)
  {
    bindex[index] = bindex[(index + 1) % 3] = min;
    bindex[(index + 2) % 3] = max;
  }
  // edge DOFs
  else
  {
    index -= 3;
    int dim = index / (order - 1);
    int offset = (index - dim * (order - 1));
    bindex[(dim + 1) % 3] = min;
    bindex[(dim + 2) % 3] = (max - 1) - offset;
    bindex[dim] = (min + 1) + offset;
  }
}


// Construct the point set in the quad element
// 1. build equispaced points with index tuple (i,j)
// 2. map index tuple to DOF index and LocalKey
template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,2>::buildQuad(std::size_t nPoints, int order, Points& points) const
{
  points.resize(nPoints);

  std::array<int,2> orders{order, order};
  std::array<Vec,4> nodes{Vec{0., 0.}, Vec{1., 0.}, Vec{1., 1.}, Vec{0., 1.}};

  for (int n = 0; n <= orders[1]; ++n) {
    for (int m = 0; m <= orders[0]; ++m) {
      // int idx = pointIndexFromIJ(m,n,orders);

      const K r = K(m) / orders[0];
      const K s = K(n) / orders[1];
      Vec p = K(1.0-r) * (nodes[3] * s + nodes[0] * K(1.0 - s)) +
              r *        (nodes[2] * s + nodes[1] * K(1.0 - s));

      auto [idx,key] = calcQuadKey(m,n,orders);
      points[idx] = LP{p, key};
      // points[idx] = LP{p, calcQuadKey(n,m,orders)};
    }
  }
}


// Obtain the VTK DOF index of the node (i,j) in the quad element
// and construct a LocalKey
template <class K>
std::pair<int,typename LagrangePointSetBuilder<K,2>::Key>
LagrangePointSetBuilder<K,2>::calcQuadKey (int i, int j, std::array<int,2> order)
{
  const bool ibdy = (i == 0 || i == order[0]);
  const bool jbdy = (j == 0 || j == order[1]);
  const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0); // How many boundaries do we lie on at once?

  int dof = 0;
  unsigned int entityIndex = 0;
  unsigned int index = 0;

  if (nbdy == 2) // Vertex DOF
  {
    dof = (i ? (j ? 2 : 1) : (j ? 3 : 0));
    entityIndex = (j ? (i ? 3 : 2) : (i ? 1 : 0));
    return std::make_pair(dof,Key{entityIndex, dim, 0});
  }

  int offset = 4;
  if (nbdy == 1) // Edge DOF
  {
    if (!ibdy) {
      dof = (i - 1) + (j ? order[0]-1 + order[1]-1 : 0) + offset;
      entityIndex = j ? 3 : 2;
      index = i-1;
    }
    else if (!jbdy) {
      dof = (j - 1) + (i ? order[0]-1 : 2 * (order[0]-1) + order[1]-1) + offset;
      entityIndex = i ? 1 : 0;
      index = j-1;
    }
    return std::make_pair(dof, Key{entityIndex, dim-1, index});
  }

  offset += 2 * (order[0]-1 + order[1]-1);

  // nbdy == 0: Face DOF
  dof = offset + (i - 1) + (order[0]-1) * ((j - 1));
  Key innerKey = LagrangePointSetBuilder<K,2>::calcQuadKey(i-1,j-1,{order[0]-2, order[1]-2}).second;
  return std::make_pair(dof, Key{0, dim-2, innerKey.index()});
}


// Lagrange points on 3d geometries
template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,3>::operator() (GeometryType gt, unsigned int order, Points& points) const
{
  std::size_t nPoints = numLagrangePoints(gt, order);

  if (gt.isTetrahedron())
    buildTetra(nPoints, order, points);
  else if (gt.isHexahedron())
    buildHex(nPoints, order, points);
  else {
    DUNE_THROW(Dune::NotImplemented,
      "Lagrange points not yet implemented for this GeometryType.");
  }

  assert(points.size() == nPoints);
}


// Construct the point set in the tetrahedron element
// 1. construct barycentric (index) coordinates
// 2. obtains the DOF index, LocalKey and actual coordinate from barycentric index
template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,3>::buildTetra (std::size_t nPoints, int order, Points& points) const
{
  points.reserve(nPoints);

  const int nVertexDOFs = 4;
  const int nEdgeDOFs = 6 * (order-1);
  const int nFaceDOFs = 4 * (order-1)*(order-2)/2;

  static const unsigned int vertexPerm[4] = {0,1,2,3};
  static const unsigned int edgePerm[6]   = {0,2,1,3,4,5};
  static const unsigned int facePerm[4]   = {1,2,0,3};

  auto calcKey = [&](int p) -> Key
  {
    if (p < nVertexDOFs) {
      return Key{vertexPerm[p], dim, 0};
    }
    else if (p < nVertexDOFs+nEdgeDOFs) {
      unsigned int entityIndex = (p - nVertexDOFs) / (order-1);
      unsigned int index = (p - nVertexDOFs) % (order-1);
      return Key{edgePerm[entityIndex], dim-1, index};
    }
    else if (p < nVertexDOFs+nEdgeDOFs+nFaceDOFs) {
      unsigned int index = (p - (nVertexDOFs + nEdgeDOFs)) % ((order-1)*(order-2)/2);
      unsigned int entityIndex = (p - (nVertexDOFs + nEdgeDOFs)) / ((order-1)*(order-2)/2);
      return Key{facePerm[entityIndex], dim-2, index};
    }
    else {
      unsigned int index = p - (nVertexDOFs + nEdgeDOFs + nFaceDOFs);
      return Key{0, dim-3, index};
    }
  };

  std::array<int,4> bindex;

  K order_d = K(order);
  for (std::size_t p = 0; p < nPoints; ++p) {
    barycentricIndex(p, bindex, order);
    Vec point{bindex[0] / order_d, bindex[1] / order_d, bindex[2] / order_d};
    points.push_back(LP{point, calcKey(p)});
  }
}


// "Barycentric index" is a set of 4 integers, each running from 0 to
// <Order>. It is the index of a point in the tetrahedron in barycentric
// coordinates.
template <class K>
void LagrangePointSetBuilder<K,3>::barycentricIndex (int p, std::array<int,4>& bindex, int order)
{
  const int nVertexDOFs = 4;
  const int nEdgeDOFs = 6 * (order-1);

  static const int edgeVertices[6][2]   = {{0,1}, {1,2}, {2,0}, {0,3}, {1,3}, {2,3}};
  static const int linearVertices[4][4] = {{0,0,0,1}, {1,0,0,0}, {0,1,0,0}, {0,0,1,0}};
  static const int vertexMaxCoords[4]   = {3,0,1,2};
  static const int faceBCoords[4][3]    = {{0,2,3}, {2,0,1}, {2,1,3}, {1,0,3}};
  static const int faceMinCoord[4]      = {1,3,0,2};

  int max = order;
  int min = 0;

  // scope into the correct tetra
  while (p >= 2 * (order * order + 1) && p != 0 && order > 3)
  {
    p -= 2 * (order * order + 1);
    max -= 3;
    min++;
    order -= 4;
  }

  // vertex DOFs
  if (p < nVertexDOFs)
  {
    for (int coord = 0; coord < 4; ++coord)
      bindex[coord] = (coord == vertexMaxCoords[p] ? max : min);
  }
  // edge DOFs
  else if (p < nVertexDOFs+nEdgeDOFs)
  {
    int edgeId = (p - nVertexDOFs) / (order-1);
    int vertexId = (p - nVertexDOFs) % (order-1);
    for (int coord = 0; coord < 4; ++coord)
    {
      bindex[coord] = min +
        (linearVertices[edgeVertices[edgeId][0]][coord] * (max - min - 1 - vertexId) +
          linearVertices[edgeVertices[edgeId][1]][coord] * (1 + vertexId));
    }
  }
  // face DOFs
  else
  {
    int faceId = (p - (nVertexDOFs+nEdgeDOFs)) / ((order-2)*(order-1)/2);
    int vertexId = (p - (nVertexDOFs+nEdgeDOFs)) % ((order-2)*(order-1)/2);

    std::array<int,3> projectedBIndex;
    if (order == 3)
      projectedBIndex[0] = projectedBIndex[1] = projectedBIndex[2] = 0;
    else
      LagrangePointSetBuilder<K,2>::barycentricIndex(vertexId, projectedBIndex, order-3);

    for (int i = 0; i < 3; i++)
      bindex[faceBCoords[faceId][i]] = (min + 1 + projectedBIndex[i]);

    bindex[faceMinCoord[faceId]] = min;
  }
}


// Construct the point set in the hexahedral element
// 1. build equispaced points with index tuple (i,j,k)
// 2. map index tuple to DOF index and LocalKey
template <class K>
  template <class Points>
void LagrangePointSetBuilder<K,3>::buildHex (std::size_t nPoints, int order, Points& points) const
{
  points.resize(nPoints);

  std::array<int,3> orders{order, order, order};
  std::array<Vec,8> nodes{Vec{0., 0., 0.}, Vec{1., 0., 0.}, Vec{1., 1., 0.}, Vec{0., 1., 0.},
                          Vec{0., 0., 1.}, Vec{1., 0., 1.}, Vec{1., 1., 1.}, Vec{0., 1., 1.}};

  for (int p = 0; p <= orders[2]; ++p) {
    for (int n = 0; n <= orders[1]; ++n) {
      for (int m = 0; m <= orders[0]; ++m) {
        const K r = K(m) / orders[0];
        const K s = K(n) / orders[1];
        const K t = K(p) / orders[2];
        Vec point = K(1.0-r) * ((nodes[3] * K(1.0-t) + nodes[7] * t) * s + (nodes[0] * K(1.0-t) + nodes[4] * t) * K(1.0-s)) +
                    r *        ((nodes[2] * K(1.0-t) + nodes[6] * t) * s + (nodes[1] * K(1.0-t) + nodes[5] * t) * K(1.0-s));

        auto [idx,key] = calcHexKey(m,n,p,orders);
        points[idx] = LP{point, key};
      }
    }
  }
}


// Obtain the VTK DOF index of the node (i,j,k) in the hexahedral element
template <class K>
std::pair<int,typename LagrangePointSetBuilder<K,3>::Key>
LagrangePointSetBuilder<K,3>::calcHexKey (int i, int j, int k, std::array<int,3> order)
{
  const bool ibdy = (i == 0 || i == order[0]);
  const bool jbdy = (j == 0 || j == order[1]);
  const bool kbdy = (k == 0 || k == order[2]);
  const int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0); // How many boundaries do we lie on at once?

  int dof = 0;
  unsigned int entityIndex = 0;
  unsigned int index = 0;

  if (nbdy == 3) // Vertex DOF
  {
    dof = (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
    entityIndex = (i ? 1 : 0) + (j ? 2 : 0) + (k ? 4 : 0);
    return std::make_pair(dof, Key{entityIndex, dim, 0});
  }

  int offset = 8;
  if (nbdy == 2) // Edge DOF
  {
    entityIndex = (k ? 8 : 4);
    if (!ibdy)
    { // On i axis
      dof = (i - 1) + (j ? order[0]-1 + order[1]-1 : 0) + (k ? 2 * (order[0]-1 + order[1]-1) : 0) + offset;
      index = i;
      entityIndex += (i ? 1 : 0);
    }
    else if (!jbdy)
    { // On j axis
      dof = (j - 1) + (i ? order[0]-1 : 2 * (order[0]-1) + order[1]-1) + (k ? 2 * (order[0]-1 + order[1]-1) : 0) + offset;
      index = j;
      entityIndex += (j ? 3 : 2);
    }
    else
    { // !kbdy, On k axis
      offset += 4 * (order[0]-1) + 4 * (order[1]-1);
      dof = (k - 1) + (order[2]-1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
      index = k;
      entityIndex = (i ? 1 : 0) + (j ? 2 : 0);
    }
    return std::make_pair(dof, Key{entityIndex, dim-1, index});
  }

  offset += 4 * (order[0]-1 + order[1]-1 + order[2]-1);
  if (nbdy == 1) // Face DOF
  {
    Key faceKey;
    if (ibdy) // On i-normal face
    {
      dof = (j - 1) + ((order[1]-1) * (k - 1)) + (i ? (order[1]-1) * (order[2]-1) : 0) + offset;
      entityIndex = (i ? 1 : 0);
      faceKey = LagrangePointSetBuilder<K,2>::calcQuadKey(j-1,k-1,{order[1]-2, order[2]-2}).second;
    }
    else {
      offset += 2 * (order[1] - 1) * (order[2] - 1);
      if (jbdy) // On j-normal face
      {
        dof = (i - 1) + ((order[0]-1) * (k - 1)) + (j ? (order[2]-1) * (order[0]-1) : 0) + offset;
        entityIndex = (j ? 3 : 2);
        faceKey = LagrangePointSetBuilder<K,2>::calcQuadKey(i-1,k-1,{order[0]-2, order[2]-2}).second;
      }
      else
      { // kbdy, On k-normal face
        offset += 2 * (order[2]-1) * (order[0]-1);
        dof = (i - 1) + ((order[0]-1) * (j - 1)) + (k ? (order[0]-1) * (order[1]-1) : 0) + offset;
        entityIndex = (k ? 5 : 4);
        faceKey = LagrangePointSetBuilder<K,2>::calcQuadKey(i-1,j-1,{order[0]-2, order[1]-2}).second;
      }
    }
    return std::make_pair(dof, Key{entityIndex, dim-2, faceKey.index()});
  }

  // nbdy == 0: Body DOF
  offset += 2 * ((order[1]-1) * (order[2]-1) + (order[2]-1) * (order[0]-1) + (order[0]-1) * (order[1]-1));
  dof = offset + (i - 1) + (order[0]-1) * ((j - 1) + (order[1]-1) * ((k - 1)));
  Key innerKey = LagrangePointSetBuilder<K,3>::calcHexKey(i-1,j-1,k-1,{order[0]-2, order[1]-2, order[2]-2}).second;
  return std::make_pair(dof, Key{0, dim-3, innerKey.index()});
}

}}} // end namespace Dune::Vtk::Impl
