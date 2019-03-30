using System;
using System.Collections.Generic;
using System.Text;
using ISAAR.MSolve.Discretization.Integration.Points;
using ISAAR.MSolve.Discretization.Integration.Quadratures;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.FEM.Interpolation;
using ISAAR.MSolve.FEM.Interpolation.GaussPointExtrapolation;
using ISAAR.MSolve.Geometry.Shapes;
using ISAAR.MSolve.Materials;

//TODO: Materials should be passed in the constructor, not each method
//TODO: not sure about the quadratures for mass. At least as many Gauss points as there are nodes are necessary 
//      (for positive definite mass matrices), but perhaps I need more for exact integration of the polynomials
//TODO: find an triangular quadrature with 6 points for Tri6 consistent mass. Symmetric triangular quadratures are the most 
//      efficient, but that doesn't concern us for mass matrices. The pdf with symmetric triangular quadratures has some other 
//      ones, as do the course notes: https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf.
namespace ISAAR.MSolve.FEM.Elements
{
    /// <summary>
    /// Creates isoparametric continuum elements with uniform thickness. Abstracts the interpolations, integrations,
    /// extrapolations and any other strategies that differentiate the elements (e.g. Quad4 from Tri6). It is also very 
    /// convenient when the thickness and material properties are the same throughout the whole domain or a region.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ContinuumElement2DFactory
    {
        private static readonly IReadOnlyDictionary<CellType, IGaussPointExtrapolation2D> extrapolations;
        private static readonly IReadOnlyDictionary<CellType, IQuadrature2D> integrationsForStiffness;
        private static readonly IReadOnlyDictionary<CellType, IQuadrature2D> integrationsForMass;
        private static readonly IReadOnlyDictionary<CellType, IIsoparametricInterpolation2D> interpolations;

        private ElasticMaterial2D_v2 commonMaterial;
        private DynamicMaterial commonDynamicProperties;
        private double commonThickness;

        static ContinuumElement2DFactory()
        {
            // Mass integrations require as many Gauss points as there are nodes, in order for the consistent mass matrix to be
            // of full rank (and symmetric positive definite)

            // Collections' declarations
            var interpolations = new Dictionary<CellType, IIsoparametricInterpolation2D>();
            var integrationsForStiffness = new Dictionary<CellType, IQuadrature2D>();
            var integrationsForMass = new Dictionary<CellType, IQuadrature2D>();
            var extrapolations = new Dictionary<CellType, IGaussPointExtrapolation2D>();

            // Quad4
            interpolations.Add(CellType.Quad4, InterpolationQuad4.UniqueInstance);
            integrationsForStiffness.Add(CellType.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
            integrationsForMass.Add(CellType.Quad4, GaussLegendre2D.GetQuadratureWithOrder(2, 2));
            extrapolations.Add(CellType.Quad4, ExtrapolationGaussLegendre2x2.UniqueInstance);

            // Quad8
            interpolations.Add(CellType.Quad8, InterpolationQuad8.UniqueInstance);
            integrationsForStiffness.Add(CellType.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            integrationsForMass.Add(CellType.Quad8, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            extrapolations.Add(CellType.Quad8, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Quad9
            interpolations.Add(CellType.Quad9, InterpolationQuad9.UniqueInstance);
            integrationsForStiffness.Add(CellType.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            integrationsForMass.Add(CellType.Quad9, GaussLegendre2D.GetQuadratureWithOrder(3, 3));
            extrapolations.Add(CellType.Quad9, ExtrapolationGaussLegendre3x3.UniqueInstance);

            // Tri3
            interpolations.Add(CellType.Tri3, InterpolationTri3.UniqueInstance);
            integrationsForStiffness.Add(CellType.Tri3, TriangleQuadratureSymmetricGaussian.Order1Point1);
            integrationsForMass.Add(CellType.Tri3, TriangleQuadratureSymmetricGaussian.Order2Points3);
            extrapolations.Add(CellType.Tri3, ExtrapolationGaussTriangular1Point.UniqueInstance);

            // Tri 6
            interpolations.Add(CellType.Tri6, InterpolationTri6.UniqueInstance);
            // see https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch24.d/IFEM.Ch24.pdf, p. 24-13, paragraph "options"
            integrationsForStiffness.Add(CellType.Tri6, TriangleQuadratureSymmetricGaussian.Order2Points3);
            integrationsForMass.Add(CellType.Tri6, TriangleQuadratureSymmetricGaussian.Order4Points6); 
            extrapolations.Add(CellType.Tri6, ExtrapolationGaussTriangular3Points.UniqueInstance);

            // Static field assignments
            ContinuumElement2DFactory.interpolations = interpolations;
            ContinuumElement2DFactory.integrationsForStiffness = integrationsForStiffness;
            ContinuumElement2DFactory.integrationsForMass = integrationsForMass;
            ContinuumElement2DFactory.extrapolations = extrapolations;
        }

        public ContinuumElement2DFactory(double commonThickness, ElasticMaterial2D_v2 commonMaterial, 
            DynamicMaterial commonDynamicProperties)
        {
            this.commonThickness = commonThickness;
            this.commonMaterial = commonMaterial;
            this.commonDynamicProperties = commonDynamicProperties;
        }

        public ContinuumElement2D CreateElement(CellType cellType, IReadOnlyList<Node_v2> nodes)
        {
            return CreateElement(cellType, nodes, commonThickness, commonMaterial, commonDynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType cellType, IReadOnlyList<Node_v2> nodes, double thickness,
            ElasticMaterial2D_v2 material, DynamicMaterial dynamicProperties)
        {
            int numGPs = integrationsForStiffness[cellType].IntegrationPoints.Count;
            var materialsAtGaussPoints = new ElasticMaterial2D_v2[numGPs];
            for (int gp = 0; gp < numGPs; ++gp) materialsAtGaussPoints[gp] = material.Clone();
            return CreateElement(cellType, nodes, thickness, materialsAtGaussPoints, dynamicProperties);
        }

        public ContinuumElement2D CreateElement(CellType cellType, IReadOnlyList<Node_v2> nodes, double thickness,
            IReadOnlyList<ElasticMaterial2D_v2> materialsAtGaussPoints, DynamicMaterial dynamicProperties)
        {
            //TODO: check if nodes - interpolation and Gauss points - materials match
            return new ContinuumElement2D(thickness, nodes, interpolations[cellType],
                integrationsForStiffness[cellType], integrationsForMass[cellType], extrapolations[cellType],
                materialsAtGaussPoints, dynamicProperties);
        }
    }
}
