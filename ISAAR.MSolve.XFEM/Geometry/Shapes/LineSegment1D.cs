using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Geometry.Shapes
{
    class LineSegment1D
    {
        public double Start { get; }
        public double End { get; }

        public LineSegment1D(double start, double end)
        {
            if (start <= end) // How to handle the degenerate case start=end?
            {
                this.Start = start;
                this.End = end;
            }
            else
            {
                this.Start = end;
                this.End = start;
            }
        }

        public SegmentSegmentPosition IntesectionWith(LineSegment1D other)
        {
            // TODO: add tolerance
            if (other.Start > this.End) return SegmentSegmentPosition.Disjoint;
            if (other.End < this.Start) return SegmentSegmentPosition.Disjoint;
            if (other.Start > this.Start && other.Start < this.End) return SegmentSegmentPosition.Overlapping; // Is the second check necessary?
            if (other.End > this.Start && other.End < this.End) return SegmentSegmentPosition.Overlapping;
            if (other.Start == this.End) return SegmentSegmentPosition.CommonVertexOnly; // TODO: return the vertex itself
            else return SegmentSegmentPosition.CommonVertexOnly; // TODO: return the vertex itself
        }

        public enum SegmentSegmentPosition
        {
            Disjoint, CommonVertexOnly, Overlapping 
        }
    }
}
