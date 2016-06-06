using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using ISAAR.MSolve.PreProcessor;
using System.Globalization;
using ISAAR.MSolve.PreProcessor.Embedding;
using ISAAR.MSolve.PreProcessor.Interfaces;
using ISAAR.MSolve.PreProcessor.Materials;
using ISAAR.MSolve.PreProcessor.Elements;

namespace ISAAR.MSolve.SimpleMeshGenerator
{
    public enum NEUDataBlock
    {
        Properties = 402,
        Nodes = 403,
        Elements = 404,
        Groups = 408,
        Constraints = 506,
        Loads = 507,
        Materials = 601
    }

    public enum NEUTopology
    {
        Line2 = 0,
        Line3,
        Tri3,
        Tri6,
        Quad4,
        Quad8,
        Tetra4,
        Wedge6,
        Brick8,
        Point,
        Tetra10,
        Wedge15,
        Brick20,
        RigidList,
        RigidList2,
        MultiList,
        Contact,
        Weld
    }

    public enum NEUElementType
    {
        Rod = 1,
        Bar,
        Tube,
        Link,
        Beam,
        Spring,
        DOFSpring,
        CurveBeam,
        Gap,
        PlotOnly,
        Shear1,
        Shear2,
        Membrane1,
        Membrane2,
        Bending1,
        Bending2,
        Plate1,
        Plate2,
        PlaneStrain1,
        PlaneStrain2,
        Laminate1,
        Laminate2,
        Axisymmetric1,
        Axisymmetric2,
        Solid1,
        Solid2,
        Mass,
        MassMat,
        Rigid,
        StiffMat,
        CurvedTube,
        PlotPlate,
        SlideLine,
        Contact,
        AxisymmetricShell1,
        AxisymmetricShell2,
        Beam82,
        Weld
    }

    public enum NEUMaterialType
    {
        Iso = 0,
        Ortho2D,
        Ortho3D,
        Aniso2D,
        Aniso3D,
        Hyperelastic,
        General,
        Fluid
    }

    public enum NEULoadType
    {
        NodalForce = 1,
        NodalMoment,
        NodalPrescribedDisplacement,
        NodalPrescribedRotation
    }

    public enum NEUEntityListType
    {
        CSys = 0,
        Point,
        Curve,
        Surface,
        Volume,
        Text,
        NotUsed,
        Node,
        Element,
        Material,
        Property,
        NodalLoad,
        ElementLoad,
        Constraint,
        ConstraintEquations,
        PointLoads,
        CurveLoads,
        SurfaceLoads,
        PointConstraints,
        CurveConstraints,
        SurfaceConstraints,
        Solids,
        ConnectionRegion,
        Connection,
        ConnectionProperty,
        Layup,
        RegionLoad
    }

    interface INEUDataBlock
    {
        void ExtractData(NEUModel model, IList<string> lines);
    }

    interface INEUProperty
    {
        int MaterialID { get; set; }
        INEUProperty BuildFromBlock(List<string> lines);
    }

    class NEUGenericProperty : INEUProperty
    {
        public int MaterialID { get; set; }

        public INEUProperty BuildFromBlock(List<string> lines)
        {
            return new NEUGenericProperty();
        }
    }

    class NEUBeamProperty : INEUProperty
    {
        public int MaterialID { get; set; }
        public double Area { get; set; }
        public double I1 { get; set; }
        public double I2 { get; set; }
        public double I12 { get; set; }
        public double J { get; set; }

        public INEUProperty BuildFromBlock(List<string> lines)
        {
            var line = lines[0].Split(',');
            return new NEUBeamProperty()
            {
                Area = Double.Parse(line[0], CultureInfo.InvariantCulture),
                I1 = Double.Parse(line[1], CultureInfo.InvariantCulture),
                I2 = Double.Parse(line[2], CultureInfo.InvariantCulture),
                I12 = Double.Parse(line[3], CultureInfo.InvariantCulture),
                J = Double.Parse(line[4], CultureInfo.InvariantCulture)
            };
        }
    }

    interface INEUMaterial
    {
        INEUMaterial BuildFromBlock(List<string> lines);
        IFiniteElementMaterial GetMaterial(NEUTopology topology);
    }

    class NEUGeneralMaterial : INEUMaterial
    {
        public INEUMaterial BuildFromBlock(List<string> lines)
        {
            return new NEUGeneralMaterial();
        }

        public IFiniteElementMaterial GetMaterial(NEUTopology topology)
        {
            throw new NotImplementedException();
        }
    }

    class NEUIsoMaterial : INEUMaterial
    {
        public double[] Es { get; set; }
        public double[] Nus { get; set; }

        public INEUMaterial BuildFromBlock(List<string> lines)
        {
            var line = lines[0].Split(',');
            return new NEUIsoMaterial()
            {
                Es = new double[] { Double.Parse(line[0], CultureInfo.InvariantCulture), Double.Parse(line[1], CultureInfo.InvariantCulture), Double.Parse(line[2], CultureInfo.InvariantCulture) },
                Nus = new double[] { Double.Parse(line[6], CultureInfo.InvariantCulture), Double.Parse(line[7], CultureInfo.InvariantCulture), Double.Parse(line[8], CultureInfo.InvariantCulture) }
            };
        }

        public IFiniteElementMaterial GetMaterial(NEUTopology topology)
        {
            switch (topology)
            {
                case NEUTopology.Brick8:
                case NEUTopology.Brick20:
                case NEUTopology.Tetra4:
                case NEUTopology.Tetra10:
                case NEUTopology.Wedge6:
                case NEUTopology.Wedge15:
                case NEUTopology.Line2:
                case NEUTopology.Line3:
                    return new ElasticMaterial3D()
                    {
                        YoungModulus = Es[0],
                        PoissonRatio = Nus[0]
                    };
                default:
                    return new ElasticMaterial()
                    {
                        YoungModulus = Es[0],
                        PoissonRatio = Nus[0]
                    };
            }
        }
    }

    class NEUModel
    {
        private readonly IDictionary<int, Node> nodes = new Dictionary<int, Node>();
        private readonly IDictionary<int, Tuple<IList<int>, int, NEUElementType, NEUTopology>> elements = new Dictionary<int, Tuple<IList<int>, int, NEUElementType, NEUTopology>>();
        private readonly IDictionary<int, INEUProperty> properties = new Dictionary<int, INEUProperty>();
        private readonly IDictionary<int, INEUMaterial> materials = new Dictionary<int, INEUMaterial>();
        private readonly IList<Tuple<int, DOFType>> constraints = new List<Tuple<int, DOFType>>();
        private readonly IList<Tuple<int, NEULoadType, int, DOFType, double>> loads = new List<Tuple<int, NEULoadType, int, DOFType, double>>();
        private readonly IDictionary<int, Tuple<string, List<int>>> groups = new Dictionary<int, Tuple<string, List<int>>>();

        public IDictionary<int, Node> Nodes { get { return nodes; } }
        public IDictionary<int, Tuple<IList<int>, int, NEUElementType, NEUTopology>> Elements { get { return elements; } }
        public IDictionary<int, INEUProperty> Properties { get { return properties; } }
        public IDictionary<int, INEUMaterial> Materials { get { return materials; } }
        public IList<Tuple<int, DOFType>> Constraints { get { return constraints; } }
        public IList<Tuple<int, NEULoadType, int, DOFType, double>> Loads { get { return loads; } }
        public IDictionary<int, Tuple<string, List<int>>> Groups { get { return groups; } }
    }

    class NEUConstraintsDataBlock : INEUDataBlock
    {
        private static readonly DOFType[] constraints = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int id = Int32.Parse(lines[2].Split(',')[0]);

            int index = 4;
            int nodeId = 0;
            while (nodeId != -1 && index < lines.Count - 1)
            {
                var nodalConstraintsLine = lines[index].Split(',');
                nodeId = Int32.Parse(nodalConstraintsLine[0]);
                if (nodeId != -1)
                {
                    for (int j = 0; j < constraints.Length; j++)
                        if (Int32.Parse(nodalConstraintsLine[3 + j]) == 1)
                            model.Constraints.Add(new Tuple<int, DOFType>(nodeId, constraints[j]));
                    index++;
                }
            }
        }
    }

    class NEUNodesDataBlock : INEUDataBlock
    {
        private static readonly DOFType[] constraints = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            for (int i = 2; i < lines.Count - 1; i++)
            {
                var separatedLine = lines[i].Split(',');
                int id = Int32.Parse(separatedLine[0]);
                var node = new Node() 
                { 
                    ID = id, 
                    X = Double.Parse(separatedLine[11], CultureInfo.InvariantCulture),
                    Y = Double.Parse(separatedLine[12], CultureInfo.InvariantCulture),
                    Z = Double.Parse(separatedLine[13], CultureInfo.InvariantCulture)
                };
                for (int j = 0; j < constraints.Length; j++)
                    if (Int32.Parse(separatedLine[5 + j]) == 1)
                        node.Constraints.Add(constraints[j]);

                model.Nodes.Add(id, node);
            }
        }
    }

    class NEUElementsDataBlock : INEUDataBlock
    {
        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int index = 2;
            while (index < lines.Count - 1)
            {
                var idLine = lines[index].Split(',');
                int id = Int32.Parse(idLine[0]);
                int propertyId = Int32.Parse(idLine[2]);
                NEUElementType typeId = (NEUElementType)Int32.Parse(idLine[3]);
                NEUTopology topologyId = (NEUTopology)Int32.Parse(idLine[4]);

                var element = new Tuple<IList<int>, int, NEUElementType, NEUTopology>(
                    lines[index + 1].Split(',').Take(10).Concat(lines[index + 2].Split(',').Take(10)).Select(x => Int32.Parse(x)).Where(x => x != 0).ToList<int>(), 
                    propertyId, typeId, topologyId);
                model.Elements.Add(id, element);

                int extraLines = 0;
                var listLine = lines[index + 6].Split(',').Take(16);
                int extraLists = listLine.Skip(12).Select(x => x != "0").Count(x => x);
                for (int i = 0; i < extraLists; i++)
                    while (!lines[index + 7 + extraLines].StartsWith("-1"))
                        extraLines++;
                index += extraLines + 7;
            }
        }
    }

    class NEUMaterialsDataBlock : INEUDataBlock
    {
        private static readonly IDictionary<NEUMaterialType, INEUMaterial> materialObjects = new Dictionary<NEUMaterialType, INEUMaterial>()
            { 
                { NEUMaterialType.Iso, new NEUIsoMaterial() }
            };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int index = 2;
            while (index < lines.Count - 1)
            {
                var idLine = lines[index].Split(',');
                int id = Int32.Parse(idLine[0]);
                NEUMaterialType typeId = (NEUMaterialType)Int32.Parse(idLine[3]);

                INEUMaterial material = materialObjects.ContainsKey(typeId) ?
                    materialObjects[typeId].BuildFromBlock(lines.ToList<string>().GetRange(index + 9, 2)) :
                    new NEUGeneralMaterial();
                model.Materials.Add(id, material);
                index += 43;
            }
        }
    }

    class NEUPropertiesDataBlock : INEUDataBlock
    {
        private static readonly IDictionary<NEUElementType, INEUProperty> propertyObjects = new Dictionary<NEUElementType, INEUProperty>()
            { 
                { NEUElementType.Rod, new NEUBeamProperty() },
                { NEUElementType.Bar, new NEUBeamProperty() },
                { NEUElementType.Beam, new NEUBeamProperty() }
            };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int index = 2;
            while (index < lines.Count - 1)
            {
                var idLine = lines[index].Split(',');
                int id = Int32.Parse(idLine[0]);
                int materialId = Int32.Parse(idLine[2]);
                NEUElementType typeId = (NEUElementType)Int32.Parse(idLine[3]);
                int laminaLines = (int)Math.Ceiling(Double.Parse(lines[index + 3], CultureInfo.InvariantCulture) / 8.0);
                int valueLines = (int)Math.Ceiling(Double.Parse(lines[index + 4 + laminaLines], CultureInfo.InvariantCulture) / 5.0);
                int outLineLines = Int32.Parse(lines[index + 5 + laminaLines + valueLines].Split(',')[0]);
                int outLine2Lines = Int32.Parse(lines[index + 6 + laminaLines + valueLines + outLineLines].Split(',')[0]);

                INEUProperty property = propertyObjects.ContainsKey(typeId) ? 
                    propertyObjects[typeId].BuildFromBlock(lines.ToList<string>().GetRange(index + 5 + laminaLines, valueLines)) :
                    new NEUGenericProperty();
                property.MaterialID = materialId;
                model.Properties.Add(id, property);

                index += 7 + laminaLines + valueLines + outLineLines + outLine2Lines;
            }
        }
    }

    class NEULoadsDataBlock : INEUDataBlock
    {
        private static readonly DOFType[] directions = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int index = 2;

            var splitLine = lines[index].Split(',');
            int loadId = Int32.Parse(splitLine[0]);
            if (lines.Count == 3) return;

            index += 21;

            // Mesh based loads
            while (!lines[index].StartsWith("-1"))
            {
                splitLine = lines[index].Split(',');
                int nodeId = Int32.Parse(splitLine[0]);
                NEULoadType loadType = (NEULoadType)Int32.Parse(splitLine[1]);
                splitLine = lines[index + 1].Split(',');
                bool[] enabledLoads = new bool[] { splitLine[0] == "1", splitLine[1] == "1", splitLine[2] == "1" };
                splitLine = lines[index + 2].Split(',');
                double[] loads = new double[] { Double.Parse(splitLine[0], CultureInfo.InvariantCulture), Double.Parse(splitLine[1], CultureInfo.InvariantCulture), Double.Parse(splitLine[2], CultureInfo.InvariantCulture) };
                for (int i = 0; i < 3; i++)
                    if (enabledLoads[i])
                        model.Loads.Add(new Tuple<int, NEULoadType, int, DOFType, double>(loadId, loadType, nodeId, directions[i], loads[i]));

                index += 7;
            }
            index++;
            // Geometry based loads
            while (!lines[index].StartsWith("-1"))
            {
                index += 19;
            }
            index++;
            // Nodal temp based loads
            while (!lines[index].StartsWith("-1"))
            {
                index += 1;
            }
            index++;
            // Element temp based loads
            while (!lines[index].StartsWith("-1"))
            {
                index += 1;
            }
            index++;
            // For version 9.3+
            //// Bolt preload loads
            //while (!lines[index].StartsWith("-1"))
            //{
            //    index += 1;
            //}
            //index++;
            //// Load definition loads
            //while (!lines[index].StartsWith("-1"))
            //{
            //    index += 1;
            //}
            //index++;
        }
    }

    class NEUGroupsDataBlock : INEUDataBlock
    {
        //private static readonly DOFType[] directions = new DOFType[] { DOFType.X, DOFType.Y, DOFType.Z, DOFType.RotX, DOFType.RotY, DOFType.RotZ };

        public void ExtractData(NEUModel model, IList<string> lines)
        {
            int index = 2;

            while (index < lines.Count - 1)
            {
                var splitLine = lines[index].Split(',');
                int groupId = Int32.Parse(splitLine[0]);
                bool needsEvaluation = splitLine[1] != "0";
                string groupName = lines[index + 1];
                if (needsEvaluation && (groupName.ToUpper().StartsWith(NEUMeshGenerator.hostGroupPrefix) || groupName.ToUpper().StartsWith(NEUMeshGenerator.embeddedGroupPrefix)))
                    throw new ArgumentException(String.Format("Group {0} needs evaluation. Run FEMAP to evaluate group, save and reload.", groupId));

                index += 24;

                // Rule types
                while (!lines[index].StartsWith("-1"))
                {
                    index++;
                    // Entries in rule type
                    while (!lines[index].StartsWith("-1"))
                    {
                        index++;
                    }
                    index++;
                }
                index++;
                // Skip line for max number of entity lists
                index++;

                while (!lines[index].StartsWith("-1"))
                {
                    splitLine = lines[index].Split(',');
                    NEUEntityListType listType = (NEUEntityListType)Int32.Parse(splitLine[0]);
                    if (listType != NEUEntityListType.Element)
                        throw new ArgumentException("Entity list type not elemental.");

                    index++;
                    List<int> ids = new List<int>();
                    while (!lines[index].StartsWith("-1"))
                    {
                        splitLine = lines[index].Split(',');
                        ids.Add(Int32.Parse(splitLine[0]));
                        index++;
                    }
                    index++;

                    model.Groups.Add(groupId, new Tuple<string, List<int>>(groupName, ids));
                }
                index++;
            }
        }
    }

    public static class NEUMeshGenerator
    {
        public const string hostGroupPrefix = "HOST_";
        public const string embeddedGroupPrefix = "EMBEDDED_";
        private const string BlockMarker = "   -1";
        private static readonly IDictionary<NEUDataBlock, INEUDataBlock> dataBlockObjects = new Dictionary<NEUDataBlock, INEUDataBlock>()
            { 
                { NEUDataBlock.Nodes, new NEUNodesDataBlock() },
                { NEUDataBlock.Elements, new NEUElementsDataBlock() },
                { NEUDataBlock.Materials, new NEUMaterialsDataBlock() },
                { NEUDataBlock.Properties, new NEUPropertiesDataBlock() },
                { NEUDataBlock.Constraints, new NEUConstraintsDataBlock() },
                { NEUDataBlock.Loads, new NEULoadsDataBlock() },
                { NEUDataBlock.Groups, new NEUGroupsDataBlock() }
            };

        public static Model ReadFile(string fileName, bool includeRotationsWhenEmbedding = false)
        {
            var model = new NEUModel();
            var fileLines = File.ReadAllLines(fileName);
            var block = new List<string>();
            IEnumerable<string> lines = null;
            lines = GetNextBlock(fileLines, block);
            while (block.Count > 0)
            {
                var blockID = GetNEUDataBlockID(block);
                if (dataBlockObjects.ContainsKey(blockID))
                    dataBlockObjects[blockID].ExtractData(model, block);

                lines = GetNextBlock(lines, block);
            }

            return GetModelFromNEUModel(model, includeRotationsWhenEmbedding);
        }

        private static void GetNEUVersion()
        {
        }

        private static Model GetModelFromNEUModel(NEUModel neuModel, bool includeRotationsWhenEmbedding)
        {
            var model = new Model();
            foreach (var node in neuModel.Nodes)
                model.NodesDictionary.Add(node.Key, node.Value);
            foreach (var constraint in neuModel.Constraints)
                model.NodesDictionary[constraint.Item1].Constraints.Add(constraint.Item2);

            foreach (var neuElement in neuModel.Elements)
            {
                var neuNodes = neuElement.Value.Item1;
                int propertyId = neuElement.Value.Item2;
                NEUElementType neuElementType = neuElement.Value.Item3;
                NEUTopology neuTopology = neuElement.Value.Item4;
                var property = neuModel.Properties[propertyId];
                var material = neuModel.Materials[property.MaterialID].GetMaterial(neuTopology);
                IFiniteElement elementType;
                switch (neuElementType)
                {
                    case NEUElementType.Solid1:
                        elementType = new Hexa8((IFiniteElementMaterial3D)material);
                        break;
                    case NEUElementType.Solid2:
                        elementType = new Hexa8((IFiniteElementMaterial3D)material);
                        //elementType = new Hexa20((IFiniteElementMaterial3D)material);
                        break;
                    case NEUElementType.Bar:
                    case NEUElementType.Beam:
                    case NEUElementType.Beam82:
                        var beamProperty = (NEUBeamProperty)property;
                        elementType = new Beam3D(material)
                        {
                            SectionArea = beamProperty.Area,
                            MomentOfInertiaZ = beamProperty.I1,
                            MomentOfInertiaY = beamProperty.I2,
                            MomentOfInertiaPolar = beamProperty.J
                        };
                        break;
                    default:
                        throw new ArgumentException("Unsupported element type.");
                }

                var element = new Element() { ElementType = elementType, ID = neuElement.Key };
                for (int i = 0; i < neuNodes.Count; i++)
                    element.AddNode(model.NodesDictionary[neuNodes[i]]);
                //element.AddNodes(model.NodesDictionary.Where(x =>
                //    {
                //        return neuNodes.IndexOf(x.Key) >= 0;
                //    }).Select(x => x.Value).ToList<Node>());
                model.ElementsDictionary.Add(neuElement.Key, element);
            }

            var dofTypeTranslations = new Dictionary<DOFType, DOFType>()
            {
                { DOFType.X, DOFType.RotX },
                { DOFType.Y, DOFType.RotY },
                { DOFType.Z, DOFType.RotZ }
            };

            foreach (var load in neuModel.Loads)
                switch (load.Item2)
                {
                    case NEULoadType.NodalForce:
                        model.Loads.Add(new Load() { Node = model.NodesDictionary[load.Item3], Amount = load.Item5, DOF = load.Item4 });
                        break;
                    case NEULoadType.NodalMoment:
                        model.Loads.Add(new Load() { Node = model.NodesDictionary[load.Item3], Amount = load.Item5, DOF = dofTypeTranslations[load.Item4] });
                        break;
                    default:
                        throw new NotImplementedException("Load type not implemented.");
                }

            var hostGroups = neuModel.Groups.Where(x => x.Value.Item1.ToUpper().StartsWith(NEUMeshGenerator.hostGroupPrefix)).Select(x => x.Value).OrderBy(x => x.Item1).ToArray();
            var embeddedGroups = neuModel.Groups.Where(x => x.Value.Item1.ToUpper().StartsWith(NEUMeshGenerator.embeddedGroupPrefix)).Select(x => x.Value).OrderBy(x => x.Item1).ToArray();
            var groupSize = hostGroups.Length;
            if (groupSize != embeddedGroups.Length)
                throw new ArgumentException("Number of host groups is not equalt to number of embedded groups. Check group names.");
            for (int i = 0; i < groupSize; i++)
            {
                var hostName = hostGroups[i].Item1.Substring(NEUMeshGenerator.hostGroupPrefix.Length);
                var embeddedName = embeddedGroups[i].Item1.Substring(NEUMeshGenerator.embeddedGroupPrefix.Length);
                if (hostName != embeddedName)
                    throw new ArgumentException(String.Format("Current host/embedded pair name does not match ({0} - {1}).", hostName, embeddedName));

                var hostElements = model.ElementsDictionary.Where(x => hostGroups[i].Item2.Any(id => id == x.Key)).Select(x => x.Value);
                var embeddedElements = model.ElementsDictionary.Where(x => embeddedGroups[i].Item2.Any(id => id == x.Key)).Select(x => x.Value);
                var embeddedGrouping = new EmbeddedGrouping(model, hostElements, embeddedElements, includeRotationsWhenEmbedding);
            }

            var subdomain = new Subdomain() { ID = 1 };
            foreach (var element in model.Elements)
                subdomain.ElementsDictionary.Add(element.ID, element);
            model.SubdomainsDictionary.Add(1, subdomain);

            return model;
        }

        private static NEUDataBlock GetNEUDataBlockID(IList<string> lines)
        {
            if (lines.Count < 3)
                throw new ArgumentException("Lines list should contain 3 or more rows.");
            return (NEUDataBlock)Int32.Parse(lines[1].Trim());
        }

        private static IEnumerable<string> GetNextBlock(IEnumerable<string> lines, List<string> linesRead)
        {
            //var linesRead = new List<string>();
            IEnumerable<string> remainingEnumerable = lines;
            linesRead.Clear();
            var insidePreviousBlockLength = remainingEnumerable.TakeWhile(x => x != BlockMarker).ToArray<string>().Length;
            remainingEnumerable = remainingEnumerable.Skip(insidePreviousBlockLength);
            var isInsidePreviousBlock = insidePreviousBlockLength > 0;
            if (isInsidePreviousBlock)
                remainingEnumerable = remainingEnumerable.Skip(1);

            var linesBeforeData = remainingEnumerable.TakeWhile(x => x == BlockMarker).ToArray<string>();
            remainingEnumerable = remainingEnumerable.Skip(linesBeforeData.Length);
            if (!isInsidePreviousBlock && linesBeforeData.Length == 0)
                return remainingEnumerable;

            if (linesBeforeData.Length != 1)
                throw new InvalidDataException("Read more than one consecutive '-1' blocks. NEU file might be corrupt.");
            linesRead.AddRange(linesBeforeData);

            var linesInsideBlock = remainingEnumerable.TakeWhile(x => x != BlockMarker).ToArray<string>();
            remainingEnumerable = remainingEnumerable.Skip(linesInsideBlock.Length);
            linesRead.AddRange(linesInsideBlock);
            if (linesRead.Count < 2)
                throw new InvalidDataException("Could not read data block identifier after '-1' block. NEU file might be corrupt.");

            var linesAfterData = remainingEnumerable.Take(1).ToArray<string>();
            remainingEnumerable = remainingEnumerable.Skip(1);
            if (linesAfterData.Length == 0 || linesAfterData[0] != BlockMarker)
                throw new InvalidDataException("Could not find '-1' blocks after reading data block. NEU file might be corrupt.");
            linesRead.AddRange(linesAfterData);
            return remainingEnumerable;
        }
    }
}
