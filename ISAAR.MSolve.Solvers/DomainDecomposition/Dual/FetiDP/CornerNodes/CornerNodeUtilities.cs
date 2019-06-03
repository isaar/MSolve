using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.Solvers.DomainDecomposition.Dual.FetiDP.CornerNodes
{
    public static class CornerNodeUtilities
    {
        public static INode[] FindCornersOfBrick3D(ISubdomain subdomain)
        {
            INode nodeXminYminZmin = null;
            INode nodeXmaxYminZmin = null;
            INode nodeXminYmaxZmin = null;
            INode nodeXmaxYmaxZmin = null;
            INode nodeXminYminZmax = null;
            INode nodeXmaxYminZmax = null;
            INode nodeXminYmaxZmax = null;
            INode nodeXmaxYmaxZmax = null;
            double xmin = Double.MaxValue;
            double xmax = Double.MinValue;
            double ymin = Double.MaxValue;
            double ymax = Double.MinValue;
            double zmin = Double.MaxValue;
            double zmax = Double.MinValue;

            foreach (INode node in subdomain.Nodes)
            {
                if (node.X < xmin) xmin = node.X;
                if (node.X > xmax) xmax = node.X;
                if (node.Y < ymin) ymin = node.Y;
                if (node.Y > ymax) ymax = node.Y;
                if (node.Z < zmin) zmin = node.Z;
                if (node.Z > zmax) zmax = node.Z;

                if (xmin == node.X && ymin == node.Y && zmin == node.Z) nodeXminYminZmin = node;
                if (xmax == node.X && ymin == node.Y && zmin == node.Z) nodeXmaxYminZmin = node;
                if (xmin == node.X && ymax == node.Y && zmin == node.Z) nodeXminYmaxZmin = node;
                if (xmax == node.X && ymax == node.Y && zmin == node.Z) nodeXmaxYmaxZmin = node;
                if (xmin == node.X && ymin == node.Y && zmax == node.Z) nodeXminYminZmax = node;
                if (xmax == node.X && ymin == node.Y && zmax == node.Z) nodeXmaxYminZmax = node;
                if (xmin == node.X && ymax == node.Y && zmax == node.Z) nodeXminYmaxZmax = node;
                if (xmax == node.X && ymax == node.Y && zmax == node.Z) nodeXmaxYmaxZmax = node;
            }

            return new INode[] 
            {
                nodeXminYminZmin, nodeXmaxYminZmin, nodeXminYmaxZmin, nodeXmaxYmaxZmin,
                nodeXminYminZmax, nodeXmaxYminZmax, nodeXminYmaxZmax, nodeXmaxYmaxZmax
            };
        }

        public static INode[] FindCornersOfRectangle2D(ISubdomain subdomain)
        {
            INode nodeXminYmin = null;
            INode nodeXmaxYmin = null;
            INode nodeXminYmax = null;
            INode nodeXmaxYmax = null;
            double xmin = Double.MaxValue;
            double xmax = Double.MinValue;
            double ymin = Double.MaxValue;
            double ymax = Double.MinValue;
            
            foreach (INode node in subdomain.Nodes)
            {
                if (node.X < xmin) xmin = node.X;
                if (node.X > xmax) xmax = node.X;
                if (node.Y < ymin) ymin = node.Y;
                if (node.Y > ymax) ymax = node.Y;

                if (xmin == node.X && ymin == node.Y) nodeXminYmin = node;
                if (xmax == node.X && ymin == node.Y) nodeXmaxYmin = node;
                if (xmin == node.X && ymax == node.Y) nodeXminYmax = node;
                if (xmax == node.X && ymax == node.Y) nodeXmaxYmax = node;
            }

            return new INode[]
            {
                nodeXminYmin, nodeXmaxYmin, nodeXminYmax, nodeXmaxYmax
            };
        }
    }
}
