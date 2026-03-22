import { test, expect } from '@playwright/test';

test.describe('HyFlow1D Web Dashboard UI Tests', () => {
  test('should display the 3D background canvas', async ({ page }) => {
    await page.goto('/');

    // Check if the Canvas element from Three.js is present
    const canvas = page.locator('canvas');
    await expect(canvas).toBeVisible();
  });

  test('should display the dashboard and run simulation button', async ({ page }) => {
    await page.goto('/');

    // Check for the CardTitle text to ensure the card rendered
    const title = page.locator('text=HyFlow1D Web Dashboard');
    await expect(title).toBeVisible();

    // Check for a specific paragraph text inside the card
    const description = page.locator('text=Hybrid FV-DG Simulation with 3D Background');
    await expect(description).toBeVisible();

    // Check if the run simulation button is rendered and visible
    const button = page.locator('button', { hasText: 'Run Simulation' });
    await expect(button).toBeVisible();
  });
});
