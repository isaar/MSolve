using System;
using ISAAR.MSolve.Discretization;
using ISAAR.MSolve.Geometry.Coordinates;

namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// Vertex of a finite element in a 1-dimensional space. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class Node1D : CartesianPoint1D, IDiscretePoint, IComparable<Node1D>
    {
        /// <summary>
        /// Instantiates a <see cref="Node1D"/>.
        /// </summary>
        /// <param name="id">A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node1D"/> 
        ///     from the rest. <see cref="ID"/> &gt;= 0.</param>
        /// <param name="x">The coordinate of the point along the single axis X.</param>
        public Node1D(int id, double x) : base(x)
        {
            if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
            this.ID = id;
        }

        /// <summary>
        /// A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node1D"/> from the rest. 
        /// <see cref="ID"/> &gt;= 0.
        /// </summary>
        public int ID { get; }

        public int CompareTo(Node1D other)
        {
            return this.ID - other.ID;
        }

        public override string ToString()
        {
            return $"node {ID} ({x})";
        }
    }
}
