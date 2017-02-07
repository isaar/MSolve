using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Geometry
{
    // The line divides the 2D space into a positive and a negative region. Points in the positive region have positive
    // signed distance to the line, while points in the negative region have negative signed distance. Each region also
    // has a normal vector. The normal of the positive region is defined by the following CONVENTION:
    // In a local 1D coordinate system, the START point of the line is to the left, the END point is to the right and
    // the normal vector of the positive region is pointing DOWNWARDS (into the positive region). 
    class Line2D : ICurve2D
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
            // The following area is positive in the positive region. TODO: prove it
            double triangleAreax2 = 
                (End.Y - Start.Y) * point.X - (End.X - Start.X) * point.Y + End.X * Start.Y - End.Y * Start.X;
            return triangleAreax2 / Length;
        }

        // The normal vector for the positive region.
        public Tuple<double, double> NormalVectorThrough(IPoint2D point)
        {
            double dy = End.Y - Start.Y;
            double dx = Start.X - End.X;
            double length = Math.Sqrt(dx * dx + dy * dy);
            return new Tuple<double, double>(dy / length, -dx / length);
        }
    }
}
