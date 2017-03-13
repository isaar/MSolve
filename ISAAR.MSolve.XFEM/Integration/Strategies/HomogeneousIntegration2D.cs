using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Materials;


namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    class HomogeneousIntegration2D: IIntegrationStrategy2D
    {
        public class Factory: IIntegrationStrategyFactory2D
        {
            private readonly IIntegrationRule2D integrationRule;
            private readonly IFiniteElementMaterial2D commonMaterial;

            public Factory(IIntegrationRule2D integrationRule, IFiniteElementMaterial2D commonMaterial)
            {
                this.integrationRule = integrationRule;
                this.commonMaterial = commonMaterial;
            }

            // So far this strategy does not need anything from the element itself. That might change imminently.
            public IIntegrationStrategy2D CreateStrategy(ContinuumElement2D element)
            {
                return new HomogeneousIntegration2D(integrationRule, commonMaterial);
            }
        }

        /// <summary>
        /// Might need it later when updating. TODO: If I had a MaterialView, I could store the same material object 
        /// for all elements and save up memory.
        /// TODO: If I need the material when updating, I will probably also need integration rule.
        /// </summary>
        private readonly IFiniteElementMaterial2D commonMaterial;
        private readonly IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> pointsAndMaterials;

        private HomogeneousIntegration2D(IIntegrationRule2D integrationRule, IFiniteElementMaterial2D commonMaterial)
        {
            this.commonMaterial = commonMaterial.Clone(); // The object passed might be mutated later on 

            IReadOnlyList<GaussPoint2D> points = integrationRule.GenerateIntegrationPoints();
            var pairs = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>[points.Count];
            for (int i = 0; i < points.Count; ++i)
            {
                pairs[i] = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>(points[i], commonMaterial.Clone());
            }
            pointsAndMaterials = pairs;
        }

        public IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> GetIntegrationPointsAndMaterials()
        {
            return pointsAndMaterials;
        }

        public void Update()
        {
            throw new NotImplementedException();
        }
    }
}
