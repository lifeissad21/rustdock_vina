import Foundation
import Metal

struct GridParams {
    var dimX: UInt32
    var dimY: UInt32
    var dimZ: UInt32
    var mapCount: UInt32
    var atomCount: UInt32
    var poseCount: UInt32
    var spacing: Float
    var slope: Float
    var originX: Float
    var originY: Float
    var originZ: Float
    var curlV: Float
}

struct PackedFloat4 {
    var x: Float
    var y: Float
    var z: Float
    var w: Float
}

struct SearchParams {
    var lanes: UInt32
    var steps: UInt32
    var seed: UInt32
    var localSteps: UInt32
    var torsionCount: UInt32
    var pairCount: UInt32
    var centerX: Float
    var centerY: Float
    var centerZ: Float
    var padding0: Float = 0
    var spanX: Float
    var spanY: Float
    var spanZ: Float
    var translationMutation: Float
    var rotationMutation: Float
    var temperature: Float
    var gradientStep: Float
    var padding1: Float = 0
}

struct TorsionData {
    var parent: UInt32
    var child: UInt32
    var maskLow: UInt32
    var maskHigh: UInt32
}

struct PairData {
    var a: UInt32
    var b: UInt32
    var typeA: UInt32
    var typeB: UInt32
}

struct Input: Decodable {
    let dims: [UInt32]
    let origin: [Float]
    let root: [Float]
    let spacing: Float
    let atoms: [[Float]]
    let torsions: [[UInt32]]
    let pairs: [[UInt32]]
    let lanes: UInt32
    let steps: UInt32
    let localSteps: UInt32
    let seed: UInt32
    let gyrationRadius: Float
    let expectedGridEnergy: Float
}
func fail(_ message: String) -> Never {
    FileHandle.standardError.write(Data((message + "\n").utf8)); exit(1)
}
func milliseconds(_ duration: Duration) -> Double {
    let c = duration.components
    return Double(c.seconds) * 1000 + Double(c.attoseconds) / 1e15
}
func run() throws {
    guard CommandLine.arguments.count == 2 else { fail("Expected a Rust-generated input directory") }
    let directory = URL(fileURLWithPath: CommandLine.arguments[1])
    let input = try JSONDecoder().decode(Input.self, from: Data(contentsOf: directory.appendingPathComponent("input.json")))
    guard input.dims.count == 3, input.origin.count == 3, input.root.count == 3,
          input.dims.allSatisfy({ $0 >= 2 }), input.spacing > 0,
          !input.atoms.isEmpty, input.atoms.count <= 64, input.atoms.allSatisfy({ $0.count == 4 }),
          input.torsions.count <= 8, input.torsions.allSatisfy({ $0.count == 4 && Int($0[0]) < input.atoms.count && Int($0[1]) < input.atoms.count }),
          input.pairs.allSatisfy({ $0.count == 4 && Int($0[0]) < input.atoms.count && Int($0[1]) < input.atoms.count }), input.lanes > 0, input.lanes <= 65536, input.steps > 0, input.localSteps > 0
    else { fail("Invalid Metal input dimensions or topology") }
    guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else { fail("No Metal compute device is available") }
    let mapData = try Data(contentsOf: directory.appendingPathComponent("maps.bin"), options: .mappedIfSafe)
    let mapStride = Int(input.dims[0]) * Int(input.dims[1]) * Int(input.dims[2]) * MemoryLayout<Float>.stride
    guard mapData.count > 0, mapData.count % mapStride == 0 else { fail("Invalid binary map size") }
    let mapCount = mapData.count / mapStride
    guard input.atoms.allSatisfy({ $0[3] >= 0 && Int($0[3]) < mapCount }) else { fail("Invalid atom map index") }
    let shader = try String(contentsOf: directory.appendingPathComponent("grid_score.metal"), encoding: .utf8)
    let options = MTLCompileOptions()
    options.fastMathEnabled = false
    let library = try device.makeLibrary(source: shader, options: options)
    func pipeline(_ name: String) throws -> MTLComputePipelineState {
        guard let function = library.makeFunction(name: name) else { fail("Missing Metal kernel \(name)") }
        return try device.makeComputePipelineState(function: function)
    }
    func buffer<T>(_ values: [T]) -> MTLBuffer {
        guard !values.isEmpty else { return device.makeBuffer(length: 16, options: .storageModeShared)! }
        return values.withUnsafeBytes { bytes in
            guard let result = device.makeBuffer(bytes: bytes.baseAddress!, length: bytes.count, options: .storageModeShared) else { fail("Metal buffer allocation failed") }
            return result
        }
    }
    func empty(_ count: Int) -> MTLBuffer {
        guard let result = device.makeBuffer(length: count * MemoryLayout<PackedFloat4>.stride, options: .storageModeShared) else { fail("Metal result buffer allocation failed") }
        return result
    }
    let maps = mapData.withUnsafeBytes { bytes in device.makeBuffer(bytes: bytes.baseAddress!, length: bytes.count, options: .storageModeShared) }
    guard let maps else { fail("Metal map allocation failed") }
    let atoms = input.atoms.map { PackedFloat4(x: $0[0], y: $0[1], z: $0[2], w: $0[3]) }
    let atomBuffer = buffer(atoms)
    var grid = GridParams(dimX: input.dims[0], dimY: input.dims[1], dimZ: input.dims[2], mapCount: UInt32(mapCount), atomCount: UInt32(atoms.count), poseCount: 1, spacing: input.spacing, slope: 1e6, originX: input.origin[0], originY: input.origin[1], originZ: input.origin[2], curlV: 1000)
    // Check this input's map packing and GPU interpolation against the Rust CPU score.
    let translations = buffer([PackedFloat4(x: input.root[0], y: input.root[1], z: input.root[2], w: 0)])
    let calibrationAtoms = buffer(atoms.map { PackedFloat4(x: 0, y: 0, z: 0, w: $0.w) })
    let calibration = empty(1)
    let scorePipeline = try pipeline("scorePoses")
    guard let check = queue.makeCommandBuffer(), let encoder = check.makeComputeCommandEncoder() else { fail("Cannot create Metal calibration command") }
    encoder.setComputePipelineState(scorePipeline)
    encoder.setBuffer(maps, offset: 0, index: 0); encoder.setBuffer(calibrationAtoms, offset: 0, index: 1)
    encoder.setBuffer(translations, offset: 0, index: 2); encoder.setBuffer(calibration, offset: 0, index: 3)
    encoder.setBytes(&grid, length: MemoryLayout<GridParams>.stride, index: 4)
    encoder.dispatchThreads(MTLSize(width: 1, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: 1, height: 1, depth: 1))
    encoder.endEncoding(); check.commit(); check.waitUntilCompleted()
    if let error = check.error { fail("Metal calibration failed: \(error)") }
    let score = calibration.contents().bindMemory(to: PackedFloat4.self, capacity: 1).pointee.w
    let gridError = abs(score - input.expectedGridEnergy)
    guard score.isFinite, gridError <= 0.002 * Float(atoms.count) else { fail("Metal map calibration differs from Rust by \(gridError) kcal/mol") }
    let lanes = Int(input.lanes)
    let steps = input.steps
    let poses = empty(lanes), orientations = empty(lanes), angles = empty(lanes * 2)
    let torsions = buffer(input.torsions.map { TorsionData(parent: $0[0], child: $0[1], maskLow: $0[2], maskHigh: $0[3]) })
    let pairs = buffer(input.pairs.map { PairData(a: $0[0], b: $0[1], typeA: $0[2], typeB: $0[3]) })
    grid.poseCount = UInt32(lanes); grid.curlV = 10
    let span = (0..<3).map { Float(input.dims[$0] - 1) * input.spacing }
    var search = SearchParams(lanes: UInt32(lanes), steps: UInt32(steps), seed: input.seed, localSteps: input.localSteps, torsionCount: UInt32(input.torsions.count), pairCount: UInt32(input.pairs.count), centerX: input.origin[0] + span[0]/2, centerY: input.origin[1] + span[1]/2, centerZ: input.origin[2] + span[2]/2, spanX: span[0], spanY: span[1], spanZ: span[2], translationMutation: 2, rotationMutation: input.gyrationRadius > 1e-6 ? 2/input.gyrationRadius : 0, temperature: 1.2, gradientStep: 0.3)
    let docking = try pipeline("flexibleDock")
    guard let command = queue.makeCommandBuffer(), let compute = command.makeComputeCommandEncoder() else { fail("Cannot create Metal docking command") }
    compute.setComputePipelineState(docking)
    for (i, b) in [maps, atomBuffer, poses, orientations].enumerated() { compute.setBuffer(b, offset: 0, index: i) }
    compute.setBytes(&grid, length: MemoryLayout<GridParams>.stride, index: 4)
    compute.setBytes(&search, length: MemoryLayout<SearchParams>.stride, index: 5)
    compute.setBuffer(torsions, offset: 0, index: 6); compute.setBuffer(pairs, offset: 0, index: 7); compute.setBuffer(angles, offset: 0, index: 8)
    let width = min(docking.threadExecutionWidth, docking.maxTotalThreadsPerThreadgroup)
    compute.dispatchThreads(MTLSize(width: lanes, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
    compute.endEncoding()
    let start = ContinuousClock.now
    command.commit(); command.waitUntilCompleted()
    if let error = command.error { fail("Metal docking failed: \(error)") }
    let elapsed = milliseconds(start.duration(to: .now))
    let p = poses.contents().bindMemory(to: PackedFloat4.self, capacity: lanes)
    let q = orientations.contents().bindMemory(to: PackedFloat4.self, capacity: lanes)
    let a = angles.contents().bindMemory(to: PackedFloat4.self, capacity: lanes * 2)
    var results = Data(capacity: lanes * 64)
    for i in 0..<lanes {
        var values = [p[i].x, p[i].y, p[i].z, p[i].w, q[i].x, q[i].y, q[i].z, q[i].w, a[2*i].x, a[2*i].y, a[2*i].z, a[2*i].w, a[2*i+1].x, a[2*i+1].y, a[2*i+1].z, a[2*i+1].w]
        values.withUnsafeMutableBytes { results.append(contentsOf: $0) }
    }
    try results.write(to: directory.appendingPathComponent("candidates.bin"), options: .atomic)
    let report: [String: Any] = ["device": device.name, "lanes": lanes, "steps": steps, "local_steps": input.localSteps, "search_ms": elapsed, "grid_error": gridError]
    print(String(data: try JSONSerialization.data(withJSONObject: report, options: .sortedKeys), encoding: .utf8)!)
}
do { try run() } catch { fail("Metal backend failed: \(error)") }
