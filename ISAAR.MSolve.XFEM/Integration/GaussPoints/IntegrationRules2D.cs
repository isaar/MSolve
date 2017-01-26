using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Integration.GaussPoints
{
    using static IntegrationRule1D;

    class IntegrationRule2D
    {
        public static readonly IntegrationRule2D Order1x1 = new IntegrationRule2D(Order1, Order1);
        public static readonly IntegrationRule2D Order2x2 = new IntegrationRule2D(Order2, Order2);
        public static readonly IntegrationRule2D Order3x3 = new IntegrationRule2D(Order3, Order3);
        public static readonly IntegrationRule2D Order4x4 = new IntegrationRule2D(Order4, Order4);

        public IReadOnlyList<GaussPoint2D> Points { get; }

        private IntegrationRule2D(IntegrationRule1D ruleXi, IntegrationRule1D ruleEta)
        {
            // Combine the integration rules of each axis: 
            var points2D = new List<GaussPoint2D>();
            foreach (var pointXi in ruleXi.Points)
            {
                foreach (var pointEta in ruleEta.Points)
                {
                    points2D.Add(new GaussPoint2D(pointXi.X, pointEta.X, pointXi.Weight * pointEta.Weight));
                }
            }
            this.Points = points2D;
        }
    }
}
