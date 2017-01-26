using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry
{
    class Line2D: ICurve2D
    {
        public IPoint2D Start { get; }
        public IPoint2D End { get; }

        public double Length
        {
            get
            {
                double dx = End.X - Start.X;
                double dy = End.Y - Start.Y;
                return Math.Sqrt(dx * dx + dy * dy);
            }
        }
        
        public Line2D(IPoint2D start, IPoint2D end)
        {
            this.Start = start;
            this.End = end;
        }

        public double SignedDistanceOf(IPoint2D point)
        {
            double triangleAreax2 = 
                (End.Y - Start.Y) * point.X - (End.X - Start.X) * point.Y + End.X * Start.Y - End.Y * Start.X;
            return triangleAreax2 / Length;
        }
    }
}
