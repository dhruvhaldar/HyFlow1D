import { serve } from "bun";
import { spawn } from "bun";

const PORT = 3001;

serve({
  port: PORT,
  async fetch(req) {
    const url = new URL(req.url);

    // Endpoint to run the C++ simulation
    if (url.pathname === "/api/run" && req.method === "POST") {
      try {
        console.log("Running simulation...");
        const proc = spawn(["../build/hyflow1d", "--clean"], {
          stdout: "pipe",
          stderr: "pipe",
        });

        await proc.exited;
        console.log("Simulation finished.");
        return new Response(JSON.stringify({ success: true }), {
          headers: { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*" },
        });
      } catch (err) {
        console.error(err);
        return new Response(JSON.stringify({ error: "Simulation failed" }), {
          status: 500,
          headers: { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*" },
        });
      }
    }

    // Endpoint to fetch results
    if (url.pathname === "/api/results" && req.method === "GET") {
      try {
        const file = Bun.file("../output/solution_final.csv");
        if (!(await file.exists())) {
          return new Response(JSON.stringify({ error: "Results not found" }), {
            status: 404,
            headers: { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*" },
          });
        }

        const text = await file.text();
        const lines = text.trim().split("\n");
        const headers = lines[0].split(",");
        const data = [];

        for (let i = 1; i < lines.length; i++) {
          if (!lines[i]) continue;
          if (lines[i].startsWith("#")) continue; // skip comments
          const values = lines[i].split(",");
          if (values.length >= 2) {
             data.push({
               x: parseFloat(values[0]),
               u: parseFloat(values[1])
             });
          }
        }

        return new Response(JSON.stringify(data), {
          headers: { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*" },
        });
      } catch (err) {
        console.error(err);
        return new Response(JSON.stringify({ error: "Failed to read results" }), {
          status: 500,
          headers: { "Content-Type": "application/json", "Access-Control-Allow-Origin": "*" },
        });
      }
    }

    // CORS preflight
    if (req.method === "OPTIONS") {
      return new Response(null, {
        headers: {
          "Access-Control-Allow-Origin": "*",
          "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
        },
      });
    }

    return new Response("Not Found", { status: 404 });
  },
});

console.log(`Backend API listening on http://localhost:${PORT}`);
