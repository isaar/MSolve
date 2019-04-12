using System;
using ISAAR.MSolve.Discretization;

namespace ISAAR.MSolve.FEM.Entities
{
    /// <summary>
    /// Vertex of a finite element in a 2-dimensional space. Immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
using System.Text;
    //public class Node2D: CartesianPoint2D, IDiscretizationPoint, IComparable<Node2D>
    public class Node2D: Node, IComparable<Node2D> // temporary, till I create the rest 2D entities
    {
        ///// <summary>
        ///// Instantiates a <see cref="Node2D"/>.
        ///// </summary>
        ///// <param name="id">A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node2D"/> 
        /////     from the rest. <see cref="ID"/> &gt;= 0.</param>
        ///// <param name="x">The coordinate of the point along axis X.</param>
        ///// <param name="y">The coordinate of the point along axis Y.</param>
        //public Node2D(int id, double x, double y): base(x, y)
        //{
        //    if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
        //    this.ID = id;
        //}

        /// <summary>
        /// Instantiates a <see cref="Node2D"/>.
        /// </summary>
        /// <param name="id">A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node2D"/> 
        ///     from the rest. <see cref="ID"/> &gt;= 0.</param>
        /// <param name="x">The coordinate of the point along axis X.</param>
        /// <param name="y">The coordinate of the point along axis Y.</param>
        public Node2D(int id, double x, double y)
        {
            if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
            this.ID = id;
            this.X = x;
            this.Y = y;
            this.Z = 0.0;
        }

        ///// <summary>
        ///// A unique identifier <see cref="ID"/> to differentiate this instance of <see cref="Node2D"/> from the rest. 
        ///// <see cref="ID"/> &gt;= 0.
        ///// </summary>
        //public int ID { get; }

        public int CompareTo(Node2D other)
        {
            return this.ID - other.ID;
        }

        //public override string ToString()
        //{
        //    return $"node {ID} ({x}, {y})";
        //}
    }
}
