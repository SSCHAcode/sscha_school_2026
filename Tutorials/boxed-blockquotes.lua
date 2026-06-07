local function contains_image(blocks)
  local found = false
  local walker = {
    Image = function()
      found = true
      return nil
    end
  }
  for _, blk in ipairs(blocks) do
    pandoc.walk_block(blk, walker)
    if found then break end
  end
  return found
end

function BlockQuote(block)
  local first = block.content[1]
  if not first or first.t ~= 'Para' then return nil end

  local first_inline = first.content[1]
  if not first_inline or first_inline.t ~= 'Strong' then return nil end

  local title = pandoc.utils.stringify(first_inline)
  local title_lower = title:lower()

  -- Skip blockquotes containing images (float figures can't live in a tcolorbox)
  if contains_image(block.content) then return nil end

  local colback, colframe, colbacktitle, coltitle
  if title_lower:match('^exercise') then
    colback, colframe, colbacktitle, coltitle = 'green!5', 'green!50!black', 'green!50!black', 'white'
  elseif title_lower:match('^question') then
    colback, colframe, colbacktitle, coltitle = 'purple!5', 'purple!60!black', 'purple!60!black', 'white'
  elseif title_lower:match('^answer') then
  colback, colframe, colbacktitle, coltitle =
    'purple!5', 'purple!60!black', 'purple!60!black', 'white'
elseif title_lower:match('^note') then
    colback, colframe, colbacktitle, coltitle = 'blue!5', 'blue!60!black', 'blue!60!black', 'white'
  else
    return nil
  end

  -- Build inner block list, removing the bold title from the first paragraph
  local inner_blocks = {}
  for i, blk in ipairs(block.content) do
    if i == 1 then
      local rest = {}
      for j = 2, #blk.content do
        rest[#rest + 1] = blk.content[j]
      end
      if #rest > 0 then
        inner_blocks[#inner_blocks + 1] = pandoc.Para(rest)
      end
    else
      inner_blocks[#inner_blocks + 1] = blk
    end
  end

  local inner_latex = ''
  if #inner_blocks > 0 then
    inner_latex = pandoc.write(pandoc.Pandoc(inner_blocks), 'latex')
  end

  local latex = '\\begin{tcolorbox}[title={\\textbf{' .. title .. '}}, colback=' .. colback .. ', colframe=' .. colframe .. ', colbacktitle=' .. colbacktitle .. ', coltitle=' .. coltitle .. ']\n' .. inner_latex .. '\\end{tcolorbox}'

  return pandoc.RawBlock('latex', latex)
end
