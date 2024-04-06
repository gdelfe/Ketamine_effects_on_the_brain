

function darkerColor = darkenColor(baseColor, darkenFactor)
    % This function takes a base RGB color and a darken factor (between 0 and 1).
    % It returns a darker version of the base color.
    % baseColor: A three-element vector specifying an RGB color.
    % darkenFactor: A scalar between 0 and 1. Higher values make the color darker.

    % Ensure the darken factor is within the range [0, 1]
    darkenFactor = max(0, min(darkenFactor, 1));

    % Reduce each color component by the darken factor
    darkerColor = baseColor * (1 - darkenFactor);
end