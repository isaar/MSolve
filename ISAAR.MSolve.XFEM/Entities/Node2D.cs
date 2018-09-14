using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.LinearAlgebra.Vectors;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;

namespace ISAAR.MSolve.XFEM.Entities
{
    // TODO: At some point I will have to override Equals() and GetHashCode(). The hash function should use 
    // the id and the number of nodes of the model
    class Node2D : ICartesianPoint2D, IComparable<Node2D>
    {
        public Node2D(int id, double x, double y)
        {
            if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
            this.ID = id;
            this.X = x;
            this.Y = y;
        }

        public Vector2 Coordinates { get { return Vector2.Create( X, Y); } }

        public int ID { get; }

        public double X { get; }

        public double Y { get; }

        public int CompareTo(Node2D other)
        {
            return this.ID - other.ID;
        }

        public override string ToString()
        {
            return "node " + ID + "(" + X + ", " + Y + ")"; 
        }
    }
}
