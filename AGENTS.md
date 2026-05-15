# AGENTS.md

## Mission

Maintain a production-grade codebase with:
- simplicity
- readability
- strong typing
- minimal dependencies
- maintainable architecture

Prefer modifying existing code over introducing abstractions.

## Workflow

For non-trivial tasks:

1. Analyze existing architecture first
2. Produce a short implementation plan
3. Implement incrementally
4. Run tests/lint/typecheck
5. Summarize changes briefly

## Coding Rules

- Prefer explicit code over magic abstractions
- Keep functions small and composable
- Avoid premature optimization
- Avoid unnecessary dependencies
- Preserve backward compatibility unless instructed otherwise
- Match existing project style

## Testing

Before completing:
- run tests
- run linter
- run type checker
- verify changed paths manually

Never claim success without verification.

## Forbidden

- Do not rewrite unrelated files
- Do not change infra/config without need
- Do not introduce frameworks
- Do not silently break APIs
- Do not suppress errors with hacks

## Architecture

- business logic must stay outside controllers/routes
- database access isolated in repositories/services
- UI components remain presentation-focused
- avoid circular dependencies

## Output Style

- concise
- technical
- no marketing language
- explain tradeoffs when relevant
