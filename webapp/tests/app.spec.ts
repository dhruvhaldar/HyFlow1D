import { test, expect } from '@playwright/test';

test.describe('Web App UI Tests', () => {
  test('should display the 3D background canvas', async ({ page }) => {
    await page.goto('/');

    // Check if the Canvas element from Three.js is present
    const canvas = page.locator('canvas');
    await expect(canvas).toBeVisible();
  });

  test('should display the glassmorphism card', async ({ page }) => {
    await page.goto('/');

    // Check for the CardTitle text to ensure the card rendered
    const title = page.locator('text=Glassmorphism UI');
    await expect(title).toBeVisible();

    // Check for a specific paragraph text inside the card
    const description = page.locator('text=Interactive 3D Background with React Three Fiber');
    await expect(description).toBeVisible();

    // Check if the button is rendered and visible
    const button = page.locator('button', { hasText: 'Get Started' });
    await expect(button).toBeVisible();
  });
});
