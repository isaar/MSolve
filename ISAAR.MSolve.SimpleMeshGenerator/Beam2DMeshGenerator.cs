using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.PreProcessor;
using ISAAR.MSolve.PreProcessor.Elements;
using ISAAR.MSolve.PreProcessor.Materials;

namespace ISAAR.MSolve.SimpleMeshGenerator
{
    public static class Beam2DMeshGenerator
    {
        public static readonly Model model = new Model();
        private static double loadAmount = 1000;
        private static double density = 7.850;
        private static double section = 0.12;
        private static double inertia = 0.0016;
        private static double young = 2.1e10;
        private static double poisson = 0.3;

        private static void GenerateNodesConstraintsAndLoads(double width, double height, double stepX, double stepY)
        {
            int nX = (int)(width / stepX) + 1;
            int nY = (int)(height / stepY) + 1;
            if (nX > nY)
            {
                for (int i = 1; i <= nX; i++)
                {
                    for (int j = 1; j <= nY; j++)
                    {
                        Node node = new Node();
                        node.ID = (i - 1) * nY + j;
                        node.X = (i - 1) * stepX;
                        node.Y = (j - 1) * stepY;
                        model.NodesDictionary.Add(node.ID, node);
                    }
                }
                for (int j = 2; j <= nY; j++)
                {
                    Load load = new Load();
                    load.Node = model.NodesDictionary[j];
                    load.DOF = DOFType.X;
                    load.Amount = loadAmount;
                    model.Loads.Add(load);
                }
            }
            else
            {
                for (int i = 1; i <= nY; i++)
                {
                    for (int j=1; j <= nX; j++)
                    {
                        Node node = new Node();
                        node.ID = (i - 1) * nX + j;
                        node.X = (j - 1) * stepX;
                        node.Y = (i - 1) * stepY;
                        model.NodesDictionary.Add(node.ID, node);
                    }
                }
                for (int i = 1; i <= nY-1; i++)
                {
                    Load load = new Load();
                    load.Node = model.NodesDictionary[i * nX + 1];
                    load.DOF = DOFType.X;
                    load.Amount = loadAmount;
                    model.Loads.Add(load);
                }
            }

            for (int i = 1; i <= nX; i++)
            {
                model.NodesDictionary[i].Constraints.Add(DOFType.X);
                model.NodesDictionary[i].Constraints.Add(DOFType.Y);
                model.NodesDictionary[i].Constraints.Add(DOFType.Z);
                model.NodesDictionary[i].Constraints.Add(DOFType.RotX);
                model.NodesDictionary[i].Constraints.Add(DOFType.RotY);
                model.NodesDictionary[i].Constraints.Add(DOFType.RotZ);
            }
        }

        private static void GenerateElements(double width, double height, double stepX, double stepY)
        {
            ElasticMaterial elasticMaterial = new ElasticMaterial();
            elasticMaterial.YoungModulus = young;
            elasticMaterial.PoissonRatio = poisson;
            Beam2D beam = new Beam2D(elasticMaterial);
            beam.Density = density;
            beam.MomentOfInertia = inertia;
            beam.SectionArea = section;
            
            int nX = (int)(width / stepX) + 1;
            int nY = (int)(height / stepY) + 1;
            int eX = 2 * nX - 1;
            int eY = nY - 1;

            for (int i = 1; i <= eY; i++)
            {
                for (int j = 1; j <= nX; j++)
                {
                    Element element = new Element();
                    element.ElementType = beam;
                    //element.MaterialType = elasticMaterial;
                    element.ID = (i - 1) * eX + (j - 1) * 2 + 1;
                    element.AddNode(model.NodesDictionary[(i - 1) * nX + j]);
                    element.AddNode(model.NodesDictionary[i * nX + j]);
                    model.ElementsDictionary.Add(element.ID, element);
                }
                for (int j = 1; j <= nX-1; j++)
                {
                    Element element = new Element();
                    element.ElementType = beam;
                    //element.MaterialType = elasticMaterial;
                    element.ID = (i - 1) * eX + (j - 1) * 2 + 2;
                    element.AddNode(model.NodesDictionary[i * nX + j]);
                    element.AddNode(model.NodesDictionary[i * nX + j + 1]);
                    model.ElementsDictionary.Add(element.ID, element);
                }
            }
        }

        private static void GenerateSubdomains(double width, double height, double stepX, double stepY, int subX, int subY)
        {
            int nX = (int)(width / stepX) + 1;
            int nY = (int)(height / stepY) + 1;
            int eX = 2 * nX - 1;
            int eY = nY - 1;
            int sX = (eX - 1) / subX;
            int sY = eY / subY;

            for (int i = 1; i <= sX; i++)
            {
                for (int j = 1; j <= sY; j++)
                {
                    Subdomain subdomain = new Subdomain();
                    subdomain.ID = (i - 1) * sY + j;
                    for (int k = 1; k <= subY; k++)
                    {
                        for (int l = 1; l <= subX; l++)
                            subdomain.ElementsDictionary.Add((i - 1) * subX + (j - 1) * subY * eX + (k - 1) * eX + l,
                                model.ElementsDictionary[(i - 1) * subX + (j - 1) * subY * eX + (k - 1) * eX + l]);
                        if (i == sX)
                            subdomain.ElementsDictionary.Add((i - 1) * subX + (j - 1) * subY * eX + (k - 1) * eX + subX + 1,
                                model.ElementsDictionary[(i - 1) * subX + (j - 1) * subY * eX + (k - 1) * eX + subX + 1]);
                    }
                    model.SubdomainsDictionary.Add(subdomain.ID, subdomain);
                }
            }
            //for (int i = 1; i <= subX; i++)
            //{
            //    for (int j = 1; j <= sY; j++)
            //    {
            //        Subdomain subdomain = new Subdomain();
            //        subdomain.ID = (i - 1) * sY + j;
            //        for (int k = 1; k <= subY; k++)
            //        {
            //            for (int l = 1; l <= sX; l++)
            //                subdomain.ElementsDictionary.Add((i - 1) * sX + (j - 1) * subY * eX + (k - 1) * eX + l, 
            //                    model.ElementsDictionary[(i - 1) * sX + (j - 1) * subY * eX + (k - 1) * eX + l]);
            //            if (i == subX)
            //                subdomain.ElementsDictionary.Add((i - 1) * sX + (j - 1) * subY * eX + (k - 1) * eX + sX + 1,
            //                    model.ElementsDictionary[(i - 1) * sX + (j - 1) * subY * eX + (k - 1) * eX + sX + 1]);
            //        }
            //        model.SubdomainsDictionary.Add(subdomain.ID, subdomain);
            //    }
            //}
        }

        private static void GenerateClusters()
        {
            Cluster cluster = new Cluster();
            cluster.ID = 1;
            foreach (Subdomain subdomain in model.Subdomains)
                cluster.Subdomains.Add(subdomain);
            model.ClustersDictionary.Add(cluster.ID, cluster);
        }

        public static void BuildMesh(double width, double height, double stepX, double stepY)
        {
            model.Clear();
            GenerateNodesConstraintsAndLoads(width, height, stepX, stepY);
            GenerateElements(width, height, stepX, stepY);
            Subdomain subdomain = new Subdomain();
            subdomain.ID = 1;
            foreach (Element element in model.Elements)
                subdomain.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary.Add(subdomain.ID, subdomain);
            GenerateClusters();
        }

        public static void BuildMesh(double width, double height, double stepX, double stepY, int subX, int subY)
        {
            model.Clear();
            GenerateNodesConstraintsAndLoads(width, height, stepX, stepY);
            GenerateElements(width, height, stepX, stepY);
            GenerateSubdomains(width, height, stepX, stepY, subX, subY);
            GenerateClusters();
        }
    }
}
