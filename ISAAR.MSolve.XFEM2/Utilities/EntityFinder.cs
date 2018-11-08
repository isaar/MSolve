using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Entities;
using ISAAR.MSolve.XFEM.Geometry.CoordinateSystems;
using ISAAR.MSolve.XFEM.Geometry.Shapes;

namespace ISAAR.MSolve.XFEM.Utilities
{
    class EntityFinder
    {
        private readonly Model2D model;
        private readonly double tolerance; // TODO: Use a value comparer instead

        public EntityFinder(Model2D model, double tolerance = 1e-4)
        {
            this.model = model;
            this.tolerance = tolerance;
        }

        /// <summary>
        /// Throws 
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <returns></returns>
        public XNode2D FindNodeWith(double x, double y)
        {
            var result = new List<XNode2D>();
            foreach (var node in model.Nodes)
            {
                if ((Math.Abs(x - node.X) <= tolerance) && (Math.Abs(y - node.Y) <= tolerance)) result.Add(node);
            }

            if (result.Count == 1) return result[0];
            else if (result.Count == 0) throw new KeyNotFoundException("No node was found, with x = " + x + " and y = " + y);
            else throw new Exception("More than 1 node were found, with x = " + x + " and y = " + y);
        }

        /// <summary>
        /// Throws 
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        public IReadOnlyList<XNode2D> FindNodesWithX(double x)
        {
            var result = new List<XNode2D>();
            foreach (var node in model.Nodes)
            {
                if (Math.Abs(x - node.X) <= tolerance) result.Add(node);
            }

            if (result.Count == 0) throw new KeyNotFoundException("No node was found, with x = " + x);
            else return result;
        }

        /// <summary>
        /// Throws 
        /// </summary>
        /// <param name="y"></param>
        /// <returns></returns>
        public IReadOnlyList<XNode2D> FindNodesWithY(double y)
        {
            var result = new List<XNode2D>();
            foreach (var node in model.Nodes)
            {
                if (Math.Abs(y - node.Y) <= tolerance) result.Add(node);
            }

            if (result.Count == 0) throw new KeyNotFoundException("No node was found, with y = " + y );
            else return result;
        }

        public List<XContinuumElement2D> FindElementsThatContains(ICartesianPoint2D point)
        {
            var result = new List<XContinuumElement2D>();
            foreach (var element in model.Elements)
            {
                var outline = ConvexPolygon2D.CreateUnsafe(element.Nodes);
                PolygonPointPosition position = outline.FindRelativePositionOfPoint(point);
                if (position == PolygonPointPosition.Inside)
                {
                    result.Add(element);
                    break;
                }
                else if ((position == PolygonPointPosition.OnEdge) || (position == PolygonPointPosition.OnVertex))
                {
                    result.Add(element);
                }
            }

            if (result.Count == 0) throw new KeyNotFoundException("No element containing the point " + point + "was found");
            return result;
        }
    }
}
