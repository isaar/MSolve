using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.Solvers.Interfaces;
using ISAAR.MSolve.FEM.Entities;
using IEmbeddedElement = ISAAR.MSolve.FEM.Interfaces.IEmbeddedElement;

namespace ISAAR.MSolve.Logging
{
    public enum NEUOutputVector
    {
        TranslationX = 2,
        TranslationY,
        TranslationZ,
        RotationX = 6,
        RotationY,
        RotationZ
    }

    public class NEUOutputVectorProperties
    {
        private readonly int id;
        private readonly string title;
        private readonly string componentVector;
        private readonly DOFType dofType;

        public int ID { get { return id; } }
        public string Title { get { return title; } }
        public string ComponentVector { get { return componentVector; } }
        public DOFType DOFType { get { return dofType; } }

        public NEUOutputVectorProperties(int id, string title, string componentVector, DOFType dofType)
        {
            this.id = id;
            this.title = title;
            this.componentVector = componentVector;
            this.dofType = dofType;
        }
    }

    public class NEUWriter
    {
        private const int outputSet = 1;
        private const string BlockMarker = "   -1";
        private readonly Dictionary<NEUOutputVector, NEUOutputVectorProperties> vectorTitles = new Dictionary<NEUOutputVector, NEUOutputVectorProperties>()
        {
            { NEUOutputVector.TranslationX, new NEUOutputVectorProperties(2, "X Translation", "2,0,0", DOFType.X) },
            { NEUOutputVector.TranslationY, new NEUOutputVectorProperties(3, "Y Translation", "0,3,0", DOFType.Y) },
            { NEUOutputVector.TranslationZ, new NEUOutputVectorProperties(4, "Z Translation", "0,0,4", DOFType.Z) },
            { NEUOutputVector.RotationX, new NEUOutputVectorProperties(6, "X Rotation", "6,0,0", DOFType.RotX) },
            { NEUOutputVector.RotationY, new NEUOutputVectorProperties(7, "Y Rotation", "0,7,0", DOFType.RotY) },
            { NEUOutputVector.RotationZ, new NEUOutputVectorProperties(8, "Z Rotation", "0,0,8", DOFType.RotZ) }
        };

        private readonly ILinearSystem subdomain;
        private readonly Model model;

        public NEUWriter(Model model, ILinearSystem subdomain)
        {
            this.model = model;
            this.subdomain = subdomain;
        }

        private List<string> MakeOutputSet(string title)
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("   450");
            lines.Add(String.Format("{0},", outputSet));
            lines.Add(title);
            lines.Add(String.Format("{0},{1},", 0, 1));
            lines.Add("0.0,");
            lines.Add("1,");
            lines.Add("<No comments>");
            lines.Add(BlockMarker);

            return lines;
        }

        private List<string> MakeNodalDataVector(NEUOutputVectorProperties vectorProperties)
        {
            var lines = new List<string>();
            lines.Add(String.Format("{0},{1},1,", outputSet, vectorProperties.ID));
            lines.Add(vectorProperties.Title);
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add(String.Format("{0},0,0,0,0,0,0,0,", vectorProperties.ComponentVector));
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("0,1,1,");

            var embeddedNodeValues = CalculateEmbeddedNodeValues();
            foreach (var node in model.NodalDOFsDictionary)
            {
                var key = new Tuple<int, DOFType>(node.Key, vectorProperties.DOFType);
                if (embeddedNodeValues.ContainsKey(key))
                    lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, embeddedNodeValues[key]));
                else
                {
                    //if (node.Value.ContainsKey(vectorProperties.DOFType) && node.Value[vectorProperties.DOFType] > -1)
                    //    lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, subdomain.Solution[node.Value[vectorProperties.DOFType]]));
                    if (!node.Value.ContainsKey(vectorProperties.DOFType)) continue;

                    if (node.Value[vectorProperties.DOFType] > -1)
                        lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, subdomain.Solution[node.Value[vectorProperties.DOFType]]));
                }
                //else
                //{
                //    var key = new Tuple<int, DOFType>(node.Key, vectorProperties.DOFType);
                //    if (embeddedNodeValues.ContainsKey(key))
                //        lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, embeddedNodeValues[key]));
                //}
            }
            lines.Add("-1,0.,");

            return lines;
        }

        private Dictionary<Tuple<int, DOFType>, double> CalculateEmbeddedNodeValues()
        {
            var embeddedNodeValues = new Dictionary<Tuple<int, DOFType>, double>();

            foreach (var element in model.Elements)
            {
                IEmbeddedElement e = element.ElementType as IEmbeddedElement;
                if (e == null) continue;

                var superElementNodes = element.ElementType.DOFEnumerator.GetNodesForMatrixAssembly(element);
                var superElementDOFs = element.ElementType.DOFEnumerator.GetDOFTypes(element);
                var superElementVector = new double[superElementDOFs.SelectMany(x => x).Count()];

                int index = 0;
                for (int i = 0; i < superElementDOFs.Count; i++)
                    for (int j = 0; j < superElementDOFs[i].Count; j++)
                    {
                        int dof = model.NodalDOFsDictionary[superElementNodes[i].ID][superElementDOFs[i][j]];
                        if (dof > -1)
                            superElementVector[index] = subdomain.Solution[dof];
                        index++;
                    }

                index = 0;
                var elementVector = e.GetLocalDOFValues(element, superElementVector);
                var elementDOFs = element.ElementType.GetElementDOFTypes(element);
                for (int i = 0; i < elementDOFs.Count; i++)
                    for (int j = 0; j < elementDOFs[i].Count; j++)
                    {
                        var key = new Tuple<int, DOFType>(element.Nodes[i].ID, elementDOFs[i][j]);
                        //int dof = model.NodalDOFsDictionary[element.Nodes[i].ID][elementDOFs[i][j]];

                        if (!embeddedNodeValues.ContainsKey(key))
                            embeddedNodeValues.Add(key, elementVector[index]);
                        index++;
                    }

                //var embeddedDOFTypes = element.ElementType.DOFEnumerator.GetDOFTypesForDOFEnumeration(element);
                //foreach (var embeddedNode in e.EmbeddedNodes)
                //{
                //    Element hostElement = embeddedNode.EmbeddedInElement;
                //    var hostDOFTypes = hostElement.ElementType.DOFEnumerator.GetDOFTypes(hostElement);
                //    var hostValues = new double[hostDOFTypes.SelectMany(x => x).Count()];
                //    int index = 0;
                //    for (int i = 0; i < hostDOFTypes.Count; i++)
                //        for (int j = 0; j < hostDOFTypes[i].Count; j++)
                //        {
                //            int dof = model.NodalDOFsDictionary[hostElement.Nodes[i].ID][hostDOFTypes[i][j]];
                //            if (dof > -1)
                //                hostValues[index] = subdomain.Solution[dof];
                //            index++;
                //        }

                //    var embeddedValues = e.GetLocalDOFValues(hostElement, hostValues);
                //    int nodeIndex = element.Nodes.IndexOf(embeddedNode.Node);
                //    index = 0;
                //    for (int i = 0; i < nodeIndex; i++)
                //        index += embeddedDOFTypes[nodeIndex].Count;
                //    //for (int i = 0; i < embeddedDOFTypes.Count; i++)
                //    for (int j = 0; j < embeddedDOFTypes[nodeIndex].Count; j++)
                //        {
                //            var key = new Tuple<int, DOFType>(element.Nodes[nodeIndex].ID, embeddedDOFTypes[nodeIndex][j]);
                //            int dof = model.NodalDOFsDictionary[element.Nodes[nodeIndex].ID][embeddedDOFTypes[nodeIndex][j]];

                //            if (dof < 0 && !embeddedNodeValues.ContainsKey(key))
                //                embeddedNodeValues.Add(key, embeddedValues[index]);
                //            index++;
                //        }
                //}
            }

            return embeddedNodeValues;
        }

        private List<string> MakeDummyTotalNodalVectors()
        {
            var lines = new List<string>();
            lines.Add(String.Format("{0},{1},1,", outputSet, 1));
            lines.Add("Total Translation");
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add("2,3,4,0,0,0,0,0,0,0,");
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("1,1,1,");
            foreach (var node in model.NodalDOFsDictionary)
            {
                double translation = 0;
                foreach (var dof in node.Value.Where(x => x.Key == DOFType.X || x.Key == DOFType.Y || x.Key == DOFType.Z).Select(x => x.Value).Where(x => x > -1))
                    translation += Math.Pow(subdomain.Solution[dof], 2);
                lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, Math.Sqrt(translation)));
            }
            lines.Add("-1,0.,");

            lines.Add(String.Format("{0},{1},1,", outputSet, 5));
            lines.Add("Total Rotation");
            lines.Add(String.Format("{0},{1},{2},", 0, 0, 0));
            lines.Add("6,7,8,0,0,0,0,0,0,0,");
            lines.Add("0,0,0,0,0,0,0,0,0,0,");
            lines.Add("0,");
            lines.Add(String.Format("{0},{1},1,7,", 0, 0));
            lines.Add("1,1,1,");
            foreach (var node in model.NodalDOFsDictionary)
            {
                double translation = 0;
                foreach (var dof in node.Value.Where(x => x.Key == DOFType.RotX || x.Key == DOFType.RotY || x.Key == DOFType.RotZ).Select(x => x.Value).Where(x => x > -1))
                    translation += Math.Pow(subdomain.Solution[dof], 2);
                lines.Add(String.Format(CultureInfo.InvariantCulture, "{0},{1},", node.Key, Math.Sqrt(translation)));
            }
            lines.Add("-1,0.,");

            return lines;
        }

        private List<string> MakeNodalDataVectors()
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("  1051");
            lines.AddRange(MakeDummyTotalNodalVectors());
            foreach (var prop in vectorTitles)
                lines.AddRange(MakeNodalDataVector(prop.Value));
            lines.Add(BlockMarker);

            return lines;
        }

        private List<string> MakeHeader()
        {
            var lines = new List<string>();
            lines.Add(BlockMarker);
            lines.Add("   100");
            lines.Add("<NULL>");
            lines.Add("10.2,");
            lines.Add(BlockMarker);
            return lines;
        }

        public void WriteToFile(string fileName)
        {
            var lines = new List<string>();
            lines.AddRange(MakeHeader());
            lines.AddRange(MakeOutputSet("Output set 1"));
            lines.AddRange(MakeNodalDataVectors());

            File.WriteAllLines(fileName, lines);
        }
    }
}
