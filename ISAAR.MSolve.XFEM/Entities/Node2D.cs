using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Entities
{
    class Node2D : IPoint2D
    {
        public int ID
        {
            get;
        }

        public double X
        {
            get;
            protected set;
        }

        public double Y
        {
            get;
            protected set;
        }

        public Node2D(int id, double x, double y)
        {
            if (id < 0) throw new ArgumentException("The parameter id must be non negative, but was: " + id);
            this.ID = id;
            this.X = x;
            this.Y = y;
        }
    }
}
