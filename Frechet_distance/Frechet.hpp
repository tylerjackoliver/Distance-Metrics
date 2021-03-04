/*  Frechet Distance
    Copyright (C) 2021 J. Tyler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef __FRECHET_H__
#define __FRECHET_H__

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <chrono>
#include "Point.hpp"

namespace Frechet
{
    /* @brief Computes the Euclidean distance between two points. The points must be of the same dimension.
       @returns The Euclidean distance between a and b; templated-type.
       @param[in] a First point.
       @param[in] b Second point.
    */
    template <typename T>
    T distanceMetric(const std::vector<T>& a, const std::vector<T>& b)
    {
        T sum = 0.0;
        for (size_t idx = 0; idx < a.size(); ++idx) sum += (a[idx] - b[idx]) * (a[idx] - b[idx]);
        return std::sqrt(sum);
    }

    /* @brief Computes the optimized distance matrix as in Devogele, T., Esnault, M., Etienne, L., & Lardy, F. (2017).
       @returns The maximum value on the 'core diagonal' of the distance matrix.
       @param[in] l1 Vector-of-Vectors containing the first trajectory.
       @param[in] l2 Vector-of-Vectors containing the second trajectory.
       @param[inout] distanceMatrix Vector-of-Vectors containing the distance between points in a and the points in b.
    */
    template <typename T>
    T computeDistanceMatrix( std::vector<std::vector<T>>& l1, std::vector<std::vector<T>>& l2, std::vector<std::vector<T>>& distanceMatrix)
    {
        /* Counters and sentinels */
        int i = 0, j = 0;
        T diagMax = 0.0;
        /* Compute trajectory lengths, remainders and quotients */
        int n = l1.size();
        int m = l2.size();
        if (n < m) /* Need to swap the trajectories */
        {
            std::swap(n, m);
            std::swap(l1, l2);
        }
        int q = static_cast<int>(n / m);
        int r = n % m;
        /* Resize the distance matrix to be of the right size. */
        distanceMatrix.reserve(n);
        for (int idx = 0; idx < n; ++idx) distanceMatrix.push_back(std::vector<T>(m));
        /* Now construct the 'almost diagonal' */
        for (i = 0; i <= r * (q + 1); ++i)
        {
            j = i / (q + 1);
            distanceMatrix[i][j] = distanceMetric(l1[i], l2[j]);
            diagMax = std::max(distanceMatrix[i][j], diagMax);
        }
        /* Construct the end of the diagonal */
        for (i = r * (q + 1); i <= (n-1); ++i)
        {
            j = (i-r) / q;
            distanceMatrix[i][j] = distanceMetric(l1[i], l2[j]);
            diagMax = std::max(distanceMatrix[i][j], diagMax);
        }
        /* Construct the upper-right matrix */
        int jMin = 0;
        for (i = 0; i <= (n-1); ++i)
        {
            j = i;
            T d;
            while ( (j < m || d <= diagMax) && j < jMin)
            {
                d = distanceMetric(l1[i], l2[j]);
                if (d < diagMax)
                {
                    distanceMatrix[i][j] = d;
                    j++;
                }
            };
            jMin = j;
        }
        /* Lower left matrix */
        int iMin = 0;
        for (j = 0; j <= (m-1); ++j)
        {
            i = j;
            T d = 0.0;
            do 
            {
                d = distanceMetric(l1[i], l2[j]);
                if (d < diagMax)
                {
                    distanceMatrix[i][j] = d;
                    i++;
                }
            } while ( (i < n || d <= diagMax) && (i < iMin) ) ;
            iMin = i;
        }
        return diagMax;
    }


    /* @brief Compute the Frechet matrix as in Devogele, T., Esnault, M., Etienne, L., & Lardy, F. (2017). Optimized Discrete Fréchet Distance between trajectories. 
       @param[in] l1 The first trajectory to compute.
       @param[in] l2 The second trajectory to compute.
       @param[inout] freechetMatrix Matrix used to determine the Frechet distance between l1 and l2.
    */
    template <typename T>
    void computeFrechetMatrix(std::vector<std::vector<T>>& l1, std::vector<std::vector<T>>& l2, std::vector<std::vector<T>>& frechetMatrix)
    {
        int i = 0, j = 0, jMin = 0;
        /* Compute the distance matrix - it will swap l1 and l2 if required,
        * so get n and m after this in case there's been a swap. */
        std::vector<std::vector<T>> distanceMatrix;
        T diagMax = computeDistanceMatrix(l1, l2, distanceMatrix);
        int n = l1.size(), m = l2.size();
        /* Resize the Frechet Matrix */
        frechetMatrix.reserve(n);
        for (unsigned idx = 0; idx < n; ++idx) frechetMatrix.push_back(std::vector<T>(m));
        frechetMatrix[0][0] = distanceMatrix[0][0];
        for (i = 1; i <= (n-1); ++i)
        {
            j = jMin;
            while ( distanceMatrix[i][j] == 0 ) j++;
            jMin = j;
            while ( distanceMatrix[i][j] != 0 && j < m)
            {
                T minimum = numeric_limits<T>::max() - 1; // Set to 'infinity'
                if (i > 0 && j > 0 && distanceMatrix[i-1][j-1] != 0)
                {
                    minimum = frechetMatrix[i-1][j-1];
                }
                if (i > 0 && distanceMatrix[i-1][j] != 0)
                {
                    minimum = std::min(minimum, frechetMatrix[i-1][j]);
                }
                if (j > 0 && distanceMatrix[i][j-1] != 0)
                {
                    minimum = std::min(minimum, frechetMatrix[i][j-1]);
                }
                frechetMatrix[i][j] = std::max(minimum, distanceMatrix[i][j]);
                j++;
            }
        }
    }

    /* @brief Computes the Frechet distance between two trajectories.
       @param[in] l1 The first trajectory to compute.
       @param[in] l2 The second trajectory to compute.
       @returns The Frechet distance between l1 and l2.
    */
    template <typename T>
    T frechetDistance(std::vector<std::vector<T>>& l1, std::vector<std::vector<T>>& l2)
    {
        std::vector<std::vector<T>> frechetMatrix;
        computeFrechetMatrix(l1, l2, frechetMatrix);
        return frechetMatrix[ l1.size() - 1 ][ l2.size() - 1 ];
    }
};
#endif
