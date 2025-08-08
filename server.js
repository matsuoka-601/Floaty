const express = require('express');
const app = express();
const path = require('path');

const PORT = 8080;

app.use((req, res, next) => {
  res.setHeader("Cross-Origin-Opener-Policy", "same-origin");
  res.setHeader("Cross-Origin-Embedder-Policy", "require-corp");
  next();
});

app.use(express.static(path.join(__dirname)));


app.listen(PORT, () => {
  console.log(`Server running at http://localhost:${PORT}/out/`);
});
