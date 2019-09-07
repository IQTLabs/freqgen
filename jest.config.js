module.exports = {
  projects: [
    {
      displayName: 'test',
      roots: ['packages/'],
      setupFilesAfterEnv: ['jest-extended'],
    },
    {
      displayName: 'lint',
      runner: 'jest-runner-eslint',
      roots: ['packages/'],
    },
  ],
  collectCoverage: true,
  collectCoverageFrom: ['packages/**/*.js', '!**/node_modules/**'],
}
