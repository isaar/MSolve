using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities;

namespace ISAAR.MSolve.XFEM.Tests
{
    static class NodeSets
    {
        public static XNode2D[] nodeSet1 = {
            new XNode2D(0, -1.0, -1.0),
            new XNode2D(1, 1.0, -1.0),
            new XNode2D(2, 1.0, 1.0),
            new XNode2D(3, -1.0, 1.0)
        };

        public static XNode2D[] nodeSet2 = {
            new XNode2D(0, -2.0, -2.0),
            new XNode2D(1, 2.0, -2.0),
            new XNode2D(2, 2.0, 2.0),
            new XNode2D(3, -2.0, 2.0)
        };

        public static XNode2D[] nodeSet3 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 2.0, 0.0),
            new XNode2D(2, 2.0, 2.0),
            new XNode2D(3, 0.0, 2.0)
        };

        public static XNode2D[] nodeSet4 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 4.0, 0.0),
            new XNode2D(2, 4.0, 4.0),
            new XNode2D(3, 0.0, 4.0)
        };

        private static XNode2D[] nodeSet5 = {
            new XNode2D(0, 0.0, 0.0),
            new XNode2D(1, 4.0, 0.0),
            new XNode2D(2, 4.0, 3.0),
            new XNode2D(3, 0.0, 3.0)
        };

        public static XNode2D[] nodeSet6 = {
            new XNode2D(0, 0.2, 0.3),
            new XNode2D(1, 2.2, 1.5),
            new XNode2D(2, 3.0, 2.7),
            new XNode2D(3, 0.7, 2.0)
        };

        public static XNode2D[] nodeSet7 = {
            new XNode2D(0, 20.0, 0.0),
            new XNode2D(1, 40.0, 0.0),
            new XNode2D(2, 40.0, 20.0),
            new XNode2D(3, 20.0, 20.0)
        };

        public static XNode2D[] nodeSet8 = {
            new XNode2D(0, 20.0, 40.0),
            new XNode2D(1, 40.0, 40.0),
            new XNode2D(2, 40.0, 60.0),
            new XNode2D(3, 20.0, 60.0)
        };
    }
}
