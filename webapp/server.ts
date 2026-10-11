import http from "node:http";
import fs from "node:fs";
import path from "node:path";
import { spawn } from "node:child_process";
import { fileURLToPath } from "node:url";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const PORT = 3001;
const ROOT_DIR = path.resolve(__dirname, "..");

interface SimulationMetrics {
  duration?: string;
  stepsPerSec?: number;
  initialMass?: number;
  finalMass?: number;
  massRelError?: number;
  energyRetention?: number;
  troubledCells?: number;
  totalDgCells?: number;
  waveRegularity?: string;
}

function parseSimulationOutput(stdout: string): SimulationMetrics {
  const metrics: SimulationMetrics = {};

  const timeMatch = stdout.match(/Simulation complete in ([\w.]+)\s+\(([\d,]+)\s+steps\/s\)/i);
  if (timeMatch) {
    metrics.duration = timeMatch[1];
    metrics.stepsPerSec = parseInt(timeMatch[2].replace(/,/g, ""), 10);
  }

  const initMassMatch = stdout.match(/Initial Mass:\s+([\d.]+)/i);
  if (initMassMatch) metrics.initialMass = parseFloat(initMassMatch[1]);

  const finalMassMatch = stdout.match(/Final Mass:\s+([\d.]+)\s+\(Rel Error:\s+([\d.]+)%\)/i);
  if (finalMassMatch) {
    metrics.finalMass = parseFloat(finalMassMatch[1]);
    metrics.massRelError = parseFloat(finalMassMatch[2]);
  }

  const energyMatch = stdout.match(/Retention:\s+([\d.]+)%/i);
  if (energyMatch) metrics.energyRetention = parseFloat(energyMatch[1]);

  const sensorMatch = stdout.match(/Troubled DG Cells:\s+(\d+)\s*\/\s*(\d+)/i);
  if (sensorMatch) {
    metrics.troubledCells = parseInt(sensorMatch[1], 10);
    metrics.totalDgCells = parseInt(sensorMatch[2], 10);
  }

  const regularityMatch = stdout.match(/Wave Regularity:\s+(.+)/i);
  if (regularityMatch) {
    metrics.waveRegularity = regularityMatch[1].trim();
  }

  return metrics;
}

const server = http.createServer(async (req, res) => {
  // CORS Headers
  res.setHeader("Access-Control-Allow-Origin", "*");
  res.setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS");
  res.setHeader("Access-Control-Allow-Headers", "Content-Type");

  if (req.method === "OPTIONS") {
    res.writeHead(204);
    res.end();
    return;
  }

  const url = new URL(req.url || "/", `http://${req.headers.host || "localhost"}`);

  // POST /api/run: Trigger simulation with parameter configuration
  if (url.pathname === "/api/run" && req.method === "POST") {
    let body = "";
    req.on("data", (chunk) => { body += chunk; });
    req.on("end", async () => {
      try {
        let params: { equation?: string; integrator?: string; threads?: number } = {};
        if (body.trim()) {
          try {
            params = JSON.parse(body);
          } catch {
            // fallback to default params
          }
        }

        const equation = params.equation === "burgers" ? "burgers" : "advection";
        const integrator = params.integrator === "euler" ? "euler" : "ssprk3";
        const threads = params.threads && params.threads > 0 ? String(params.threads) : "4";

        const binaryPath = path.join(ROOT_DIR, "build", "hyflow1d");
        const args = ["--clean", "-e", equation, "-i", integrator, "-t", threads, "--vtk"];

        console.log(`Running simulation: ${binaryPath} ${args.join(" ")}`);

        const child = spawn(binaryPath, args, {
          cwd: ROOT_DIR,
          env: { ...process.env, NO_COLOR: "1" }
        });

        let stdout = "";
        let stderr = "";

        child.stdout.on("data", (data) => { stdout += data.toString(); });
        child.stderr.on("data", (data) => { stderr += data.toString(); });

        child.on("close", (code) => {
          if (code === 0) {
            const metrics = parseSimulationOutput(stdout);
            res.writeHead(200, { "Content-Type": "application/json" });
            res.end(JSON.stringify({ success: true, metrics, stdout }));
          } else {
            console.error("Simulation failed:", stderr || stdout);
            res.writeHead(500, { "Content-Type": "application/json" });
            res.end(JSON.stringify({ error: "Simulation failed", stderr, stdout }));
          }
        });
      } catch (err: any) {
        console.error("Server error during run:", err);
        res.writeHead(500, { "Content-Type": "application/json" });
        res.end(JSON.stringify({ error: err.message || "Internal server error" }));
      }
    });
    return;
  }

  // GET /api/results: Fetch latest solution curve
  if (url.pathname === "/api/results" && req.method === "GET") {
    try {
      const candidatePaths = [
        path.join(ROOT_DIR, "output", "solution_final.csv"),
        path.join(ROOT_DIR, "build", "output", "solution_final.csv")
      ];

      let solutionPath = "";
      for (const p of candidatePaths) {
        if (fs.existsSync(p)) {
          solutionPath = p;
          break;
        }
      }

      if (!solutionPath) {
        res.writeHead(404, { "Content-Type": "application/json" });
        res.end(JSON.stringify({ error: "Results not found. Run simulation first." }));
        return;
      }

      const text = fs.readFileSync(solutionPath, "utf-8");
      const lines = text.trim().split("\n");
      const data: Array<{ x: number; u: number }> = [];
      const meta: Record<string, string> = {};

      for (let i = 0; i < lines.length; i++) {
        const line = lines[i].trim();
        if (!line) continue;
        if (line.startsWith("#")) {
          const content = line.substring(1).trim();
          const tokens = content.split(/\s+/);
          for (const token of tokens) {
            const [k, v] = token.split("=");
            if (k && v) meta[k] = v;
          }
          continue;
        }
        const values = line.split(",");
        if (values.length >= 2 && values[0].toLowerCase() !== "x") {
          const x = parseFloat(values[0]);
          const u = parseFloat(values[1]);
          if (!isNaN(x) && !isNaN(u)) {
            data.push({ x, u });
          }
        }
      }

      // Backward compatible: return array directly if requested or rich object
      res.writeHead(200, { "Content-Type": "application/json" });
      res.end(JSON.stringify(data));
      return;
    } catch (err: any) {
      console.error("Error reading results:", err);
      res.writeHead(500, { "Content-Type": "application/json" });
      res.end(JSON.stringify({ error: "Failed to read results" }));
      return;
    }
  }

  res.writeHead(404, { "Content-Type": "text/plain" });
  res.end("Not Found");
});

server.listen(PORT, () => {
  console.log(`HyFlow1D Backend API listening on http://localhost:${PORT}`);
});
